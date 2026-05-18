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


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--n-indiv",     type=int,   default=200,  help="diploid individuals (training set)")
    ap.add_argument("--n-test",      type=int,   default=0,    help="additional held-out test individuals (0 = no test split)")
    ap.add_argument("--length",      type=float, default=5e6,  help="sequence length (bp)")
    ap.add_argument("--recomb-rate", type=float, default=1e-8, help="recombination rate per bp per gen")
    ap.add_argument("--mut-rate",    type=float, default=1e-8, help="mutation rate per bp per gen")
    ap.add_argument("--Ne",          type=float, default=1e4,  help="effective population size")
    ap.add_argument("--h2-fixed",    type=float, default=0.3,  help="target fixed-effect heritability")
    ap.add_argument("--h2-random",   type=float, default=0.3,  help="target random-effect heritability")
    ap.add_argument("--n-causal",    type=int,   default=20,   help="number of causal SNPs (fixed effects)")
    ap.add_argument("--seed",        type=int,   default=42,   help="random seed")
    ap.add_argument("--out",         type=str,   default="sim", help="output prefix")
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

    if args.h2_fixed + args.h2_random > 1.0:
        sys.exit("ERROR: h2_fixed + h2_random must be ≤ 1")
    if args.h2_fixed < 0 or args.h2_random < 0:
        sys.exit("ERROR: heritabilities must be non-negative")

    rng = np.random.default_rng(args.seed)

    # ---- 1. Simulate ARG and mutations -------------------------------------
    n_total = args.n_indiv + args.n_test
    print(f"Simulating ARG: n_train={args.n_indiv} + n_test={args.n_test} "
          f"= {n_total} individuals, L={args.length:.0f} bp ...", flush=True)
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

    p = len(dosage_list)
    if p == 0:
        sys.exit("ERROR: no polymorphic biallelic variants found")
    print(f"  {p} polymorphic variants", flush=True)

    dosage = np.stack(dosage_list)   # (p, n_indiv)
    positions = pos_list

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

    n_causal = min(args.n_causal, p)
    causal_idx = rng.choice(p, size=n_causal, replace=False)
    causal_idx.sort()

    beta = np.zeros(p)
    for j in causal_idx:
        k = rng.choice(K, p=pi_mix)
        if sigma2_mix[k] > 0:
            beta[j] = rng.normal(0.0, math.sqrt(sigma2_mix[k]))

    # Fixed-effect phenotype contribution (uses standardized genotypes)
    y_fixed = dosage_std.T @ beta           # (n_indiv,)

    var_fixed = float(np.var(y_fixed))
    if var_fixed > 0 and args.h2_fixed > 0:
        scale_f = math.sqrt(args.h2_fixed / var_fixed)
        beta   *= scale_f
        y_fixed = dosage_std.T @ beta
    elif args.h2_fixed == 0:
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

    # Draw ws² for every tree from the discrete mixture
    mix_draw = rng.choice(M, size=n_trees, p=phi_mix)
    ws2_vec  = tau2_mix[mix_draw]       # shape (n_trees,)

    print("Simulating random effects ...", flush=True)
    y_random = sim_random_effect(ts, hap_ranks, n_indiv, ws2_vec, rng)

    var_random = float(np.var(y_random))
    if var_random > 0 and args.h2_random > 0:
        scale_r   = math.sqrt(args.h2_random / var_random)
        y_random *= scale_r
        ws2_vec  *= scale_r ** 2        # scale the variance parameters accordingly
    elif args.h2_random == 0:
        y_random[:] = 0.0
        ws2_vec[:]  = 0.0

    # ---- 5. Residual noise and final phenotype ----------------------------
    sigma2 = max(0.0, 1.0 - args.h2_fixed - args.h2_random)
    y_noise = rng.normal(0.0, math.sqrt(sigma2) if sigma2 > 0 else 0.0, size=n_indiv)
    y = y_fixed + y_random + y_noise

    # Realised heritabilities (from rescaled components, before noise is added)
    h2_fixed_realised  = float(np.var(y_fixed)  / np.var(y)) if np.var(y) > 0 else 0.0
    h2_random_realised = float(np.var(y_random) / np.var(y)) if np.var(y) > 0 else 0.0

    # ---- 6. Write outputs -------------------------------------------------
    print(f"Writing outputs with prefix '{args.out}' ...", flush=True)

    n_train = args.n_indiv
    n_test  = args.n_test
    has_test = n_test > 0
    ts_train = None

    if has_test:
        # Save full ARG (needed for cross-kinship between train and test in prediction)
        ts.dump(args.out + ".full.trees")
        # Build the train tree sequence by simplifying to train sample nodes
        train_sample_nodes = []
        for i in range(n_train):
            train_sample_nodes.extend(ts.individual(i).nodes)
        ts_train = ts.simplify(samples=train_sample_nodes)
        ts_train.dump(args.out + ".trees")

        # Train: phenotypes, genotypes for first n_train individuals
        np.savetxt(args.out + ".phen", y[:n_train], fmt="%.8g")
        write_plink(args.out, dosage[:, :n_train], positions, n_train, p)

        # Test: phenotypes, genotypes for last n_test individuals
        np.savetxt(args.out + ".test.phen", y[n_train:], fmt="%.8g")
        write_plink(args.out + ".test", dosage[:, n_train:], positions, n_test, p)

        # Also save full phen for completeness
        np.savetxt(args.out + ".full.phen", y, fmt="%.8g")
        # Per-individual true components (for evaluation of prediction quality)
        with open(args.out + ".test.true_components.tsv", "w") as f:
            f.write("idx\ty_fixed\ty_random\ty_total\n")
            for i in range(n_test):
                gi = n_train + i
                f.write(f"{i}\t{y_fixed[gi]:.8g}\t{y_random[gi]:.8g}\t{y[gi]:.8g}\n")
    else:
        ts.dump(args.out + ".trees")
        np.savetxt(args.out + ".phen", y, fmt="%.8g")
        write_plink(args.out, dosage, positions, n_train, p)

    with open(args.out + ".true_beta.tsv", "w") as f:
        f.write("snp_idx\tbp_pos\tbeta\n")
        for j in causal_idx:
            f.write(f"{j}\t{positions[j]}\t{beta[j]:.8g}\n")

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

    params = {
        "n_indiv":         n_train,
        "n_test":          n_test,
        "n_total":         n_train + n_test,
        "n_trees":         (ts_train.num_trees if (has_test and ts_train is not None)
                            else n_trees),
        "n_trees_full":    n_trees,
        "n_variants":      p,
        "n_causal":        n_causal,
        "length":          args.length,
        "recomb_rate":     args.recomb_rate,
        "mut_rate":        args.mut_rate,
        "Ne":              args.Ne,
        "seed":            args.seed,
        "h2_fixed_target":  args.h2_fixed,
        "h2_random_target": args.h2_random,
        "h2_fixed_realised":  h2_fixed_realised,
        "h2_random_realised": h2_random_realised,
        "sigma2":          sigma2,
        "regime":          args.regime,
        "sigma2_mix":      list(sigma2_mix),
        "pi_mix":          list(pi_mix),
        "tau2_frac":       list(tau2_frac),
        "tau2_mix":        list(tau2_mix),
        "phi_mix":         list(phi_mix),
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
    print(f"  h2_fixed  target / realised : {args.h2_fixed:.3f} / {h2_fixed_realised:.4f}")
    print(f"  h2_random target / realised : {args.h2_random:.3f} / {h2_random_realised:.4f}")
    print(f"  sigma2           : {sigma2:.4f}")


if __name__ == "__main__":
    main()
