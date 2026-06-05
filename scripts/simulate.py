#!/usr/bin/env python3
"""
simulate.py — generate synthetic ARG + genotype + phenotype data for SARGENT testing.

Model:  y = X β + u + ε
  X β : fixed effects (BayesR discrete mixture)
   u  : random effects drawn from tree-branch Brownian motion with per-tree ws²
   ε  : i.i.d. N(0, σ²)
   V  = Σ_s ws²[s] G_s + σ² I

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

def sim_random_effect(ts, hap_ranks, n_indiv, ws2_vec, rng):
    """Return u vector (length n_indiv) for all trees.

    ws2_vec: array of ws² values, one per tree (length ts.num_trees).
    hap_ranks: shape (n_indiv, 2) — sample ranks for each individual's two haplotypes.
    """
    u = np.zeros(n_indiv)
    node_val = np.zeros(ts.num_nodes)

    for tree_idx, tree in enumerate(ts.trees()):
        ws2 = ws2_vec[tree_idx]
        if ws2 == 0.0:
            continue

        node_val[:] = 0.0
        # Preorder: root first, then children — parent value is ready before children
        for nd in tree.nodes(order="preorder"):
            parent = tree.parent(nd)
            if parent == tskit.NULL:
                node_val[nd] = 0.0
            else:
                bl = tree.branch_length(nd)
                node_val[nd] = node_val[parent] + rng.normal(0.0, math.sqrt(ws2 * bl))

        # Accumulate into u: (1/√2) · (val[hap0] + val[hap1])
        h0 = hap_ranks[:, 0]
        h1 = hap_ranks[:, 1]
        # hap_ranks are sample ranks; convert to node IDs via ts.samples()
        samples = ts.samples()
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
                   Requires `ts` (tskit.TreeSequence). Designed for Fig 5 Panel B
                   — gives maximum age-bin separability (pure 0/1 truth profile).

    Absolute scale is arbitrary — caller renormalises ws2_vec so realised var(u)
    matches h2_arg; only the relative profile matters downstream.
    """
    if n_trees <= 0:
        return np.zeros(0, dtype=float)
    if kind == "diffuse":
        return np.ones(n_trees, dtype=float)

    if kind.startswith("age_"):
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
                    choices=["panmictic", "admixed", "transfer"],
                    help="population model. 'panmictic' = single constant-size pop; "
                         "'admixed' = two diverged sources + recent admixture pulse, "
                         "giving local-ancestry heterogeneity along the genome (the "
                         "Phase 10 optimistic regime, where different regions carry "
                         "genuinely different relatedness so per-segment ws² can help); "
                         "'transfer' = two diverged pops A and B (no admixture pulse) "
                         "with training drawn from A + optional B fraction (--mix-fraction-B) "
                         "and test drawn entirely from B (Block 1.1 cross-ancestry PRS).")
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
    ap.add_argument("--n-causal",    type=int,   default=None,
                    help="number of causal SNPs (fixed effects).  If unset, defaults to "
                         "max(20, int(--causal-fraction * n_variants)).  An explicit value "
                         "takes priority over --causal-fraction.")
    ap.add_argument("--causal-fraction", type=float, default=0.01,
                    help="fraction of polymorphic variants that are causal SNPs "
                         "(default 0.01 = 1%%, matches SBayesR's canonical scenario for human "
                         "polygenic traits; typical GWAS-polygenic range is 0.01-0.10).  "
                         "Overridden by an explicit --n-causal.")
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
                         "what creates the flavour-A-vs-B gap).")
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
                             "age_recent", "age_middle", "age_ancient"],
                    help="directly build the per-tree ws² profile (Phase 10 / Block 2.6), "
                         "overriding the --regime mixture draw.  "
                         "'diffuse' = flat constant ws²; "
                         "'sparse' = floor + n_blocks evenly-spaced contiguous blocks; "
                         "'scattered' = floor + block_frac·n_trees individually scattered "
                         "active trees (n_blocks ignored). The 'scattered' option is the "
                         "fused-lasso robustness check: confirms that any per-segment "
                         "advantage isn't an artifact of the prior matching contiguous truth.")
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
    if args.demography == "transfer":
        n_train_B = int(round(args.mix_fraction_B * args.n_indiv))
        n_train_A = args.n_indiv - n_train_B

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
    else:
        ts = msprime.sim_ancestry(
            samples=n_total,
            sequence_length=args.length,
            recombination_rate=args.recomb_rate,
            population_size=args.Ne,
            random_seed=args.seed,
        )
    ts = msprime.sim_mutations(ts, rate=args.mut_rate, random_seed=args.seed + 1)
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

    p = len(dosage_list)
    if p == 0:
        sys.exit("ERROR: no polymorphic biallelic variants found")
    print(f"  {p} polymorphic variants", flush=True)

    dosage = np.stack(dosage_list)   # (p, n_indiv)
    positions = pos_list
    ages = np.array(age_list, dtype=float)  # per-variant age in generations

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
    if args.h2_unobs > 0:
        n_unobs_causal = min(args.n_unobs_causal, p)
        freq = dosage.mean(axis=1) / 2.0
        maf  = np.minimum(freq, 1.0 - freq)
        eligible = np.where((maf >= args.unobs_maf_min) &
                            (maf <= args.unobs_maf_max))[0]
        if len(eligible) < n_unobs_causal:
            sys.exit(f"ERROR: only {len(eligible)} variants in MAF "
                     f"[{args.unobs_maf_min}, {args.unobs_maf_max}], need "
                     f"{n_unobs_causal}.  Lengthen --length or widen MAF range.")
        if var_tau is not None:
            ew = var_tau[eligible].astype(float)
            unobs_idx = rng.choice(eligible, size=n_unobs_causal, replace=False,
                                   p=ew / ew.sum())
        else:
            unobs_idx = rng.choice(eligible, size=n_unobs_causal, replace=False)
        unobs_idx.sort()
        unobs_set = set(unobs_idx.tolist())
        remaining = np.array([j for j in range(p) if j not in unobs_set])
    else:
        n_unobs_causal = 0
        unobs_idx = np.zeros(0, dtype=int)
        remaining = np.arange(p)

    if args.bed_coverage < 1.0:
        bed_rng = np.random.default_rng(args.seed + 1_000_003)
        bed_mask = bed_rng.random(p) < args.bed_coverage
    else:
        bed_mask = np.ones(p, dtype=bool)

    # Phase 6b: unobserved-class variants are never in the BED, regardless of
    # coverage f.  Force-exclude them after the f-coverage draw.
    if args.h2_unobs > 0:
        bed_mask[unobs_idx] = False

    # Block 2.2: common-SNP / array-like typed panel.  Drop rare variants from
    # the BED so they are untyped (only the ARG can see them).
    if args.bed_maf_min > 0.0:
        freq_all = dosage.mean(axis=1) / 2.0
        maf_all  = np.minimum(freq_all, 1.0 - freq_all)
        bed_mask &= (maf_all >= args.bed_maf_min)

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
    if args.n_causal is not None:
        n_causal_target = min(args.n_causal, len(remaining))
    else:
        n_causal_target = min(max(20, int(args.causal_fraction * p)), len(remaining))
    n_causal = n_causal_target
    # Block 2.2 — annotation enrichment via causal DENSITY: a variant's causal
    # probability is proportional to tau_{ann(j)}, so enriched annotations carry
    # tau-fold more causal variants per base (constant effect sizes).  This makes
    # per-annotation h^2 proportional to base_frac * tau while keeping the signal
    # polygenic (many modest-effect SNPs), which the downstream estimators need.
    if var_tau is not None:
        w = var_tau[remaining].astype(float)
        causal_idx = rng.choice(remaining, size=n_causal, replace=False,
                                p=w / w.sum())
    else:
        causal_idx = rng.choice(remaining, size=n_causal, replace=False)
    causal_idx.sort()

    beta = np.zeros(p)
    for j in causal_idx:
        k = rng.choice(K, p=pi_mix)
        if sigma2_mix[k] > 0:
            beta[j] = rng.normal(0.0, math.sqrt(sigma2_mix[k]))

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

    var_fixed = float(np.var(y_fixed))
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

    var_random = float(np.var(y_random))
    if var_random > 0 and args.h2_arg > 0:
        scale_r   = math.sqrt(args.h2_arg / var_random)
        y_random *= scale_r
        ws2_vec  *= scale_r ** 2        # scale the variance parameters accordingly
    elif args.h2_arg == 0:
        y_random[:] = 0.0
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
        beta_unobs[unobs_idx] = eff
        y_unobs = dosage_std.T @ beta_unobs

        var_unobs = float(np.var(y_unobs))
        if var_unobs > 0:
            scale_u = math.sqrt(args.h2_unobs / var_unobs)
            beta_unobs *= scale_u
            y_unobs    *= scale_u

    # ---- 5. Residual noise and final phenotype ----------------------------
    sigma2 = max(0.0, 1.0 - args.h2_snp - args.h2_arg - args.h2_unobs)
    y_noise = rng.normal(0.0, math.sqrt(sigma2) if sigma2 > 0 else 0.0, size=n_indiv)
    y = y_fixed + y_random + y_unobs + y_noise

    # Realised heritabilities (from rescaled components, before noise is added)
    h2_snp_realised   = float(np.var(y_fixed)  / np.var(y)) if np.var(y) > 0 else 0.0
    h2_arg_realised   = float(np.var(y_random) / np.var(y)) if np.var(y) > 0 else 0.0
    h2_unobs_realised = float(np.var(y_unobs)  / np.var(y)) if np.var(y) > 0 else 0.0

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
    else:
        ts.dump(args.out + ".trees")
        np.savetxt(args.out + ".phen", y, fmt="%.8g")
        write_plink(args.out, bed_dosage, bed_positions, n_train, n_bed)

    # Map full-panel indices to BED-panel indices (for SNPs that are in the BED)
    bed_idx_of = np.full(p, -1, dtype=int)
    bed_idx_of[bed_mask] = np.arange(n_bed)

    with open(args.out + ".true_beta.tsv", "w") as f:
        f.write("snp_idx\tbp_pos\tbeta\tin_bed\tbed_idx\tage\n")
        for j in causal_idx:
            in_bed = 1 if bed_mask[j] else 0
            f.write(f"{j}\t{positions[j]}\t{beta[j]:.8g}\t{in_bed}\t"
                    f"{bed_idx_of[j] if in_bed else -1}\t{ages[j]:.6g}\n")

    # ws2_vec is defined on the FULL ARG's trees (or on the only ARG if no test split).
    # If there is a test split, we also record train-tree ws² indexed by the simplified
    # train tree sequence. Each train tree inherits ws² from the full-ARG tree containing
    # its midpoint.
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
        "n_causal":        n_causal,
        "n_causal_in_bed": int(bed_mask[causal_idx].sum()),
        "length":          args.length,
        "recomb_rate":     args.recomb_rate,
        "mut_rate":        args.mut_rate,
        "Ne":              args.Ne,
        "demography":      args.demography,
        "admix_split_time": args.admix_split_time if args.demography == "admixed" else None,
        "admix_time":       args.admix_time       if args.demography == "admixed" else None,
        "admix_alpha":      args.admix_alpha      if args.demography == "admixed" else None,
        "transfer_split_time": args.transfer_split_time if args.demography == "transfer" else None,
        "mix_fraction_B":     args.mix_fraction_B     if args.demography == "transfer" else None,
        "n_train_A":          n_train_A               if args.demography == "transfer" else None,
        "n_train_B":          n_train_B               if args.demography == "transfer" else None,
        "seed":            args.seed,
        "h2_snp_target":   args.h2_snp,
        "h2_arg_target":   args.h2_arg,
        "h2_unobs_target": args.h2_unobs,
        "age_tau_coef":    args.age_tau_coef,
        "age_tau_ref":     args.age_tau_ref if args.age_tau_coef != 0 else None,
        "h2_snp_realised":   h2_snp_realised,
        "h2_arg_realised":   h2_arg_realised,
        "h2_unobs_realised": h2_unobs_realised,
        "n_unobs_causal":    n_unobs_causal,
        "unobs_effect_dist":  args.unobs_effect_dist if args.h2_unobs > 0 else None,
        "unobs_effect_scale": args.unobs_effect_scale if args.h2_unobs > 0 else None,
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
