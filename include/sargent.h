/*
 * SARGENT: Sparse ARG-based ENcoding of Traits
 *
 * Copyright (c) 2026, Cold Spring Harbor Laboratory
 * All rights reserved.
 *
 * This file is part of SARGENT and is distributed under the BSD 3-Clause License.
 * See the LICENSE file in the project root for details.
 */

#pragma once

#include <tskit.h>
#include "geno.h"
#include "gall.h"

typedef struct {
    int n;              /* number of individuals */
    int p;              /* number of markers */
    int S;              /* number of segments (local trees) */
    geno_reader_t  *gr; /* genotype reader (borrowed) */
    tsk_treeseq_t  *ts; /* tree sequence (borrowed) */
    double         *y;  /* phenotype vector [n] */
    double        *beta;   /* SNP effects [p], NULL until inferred */
    double        *ws2;    /* segment weights [S], NULL until inferred */
    double h2;          /* total heritability (h2_snp + h2_arg) */
    double h2_snp;      /* variance explained by BayesR fit to typed SNPs */
    double h2_arg;      /* additional variance captured by the ARG, after the
                           SNP fit; "missing heritability recovered" */
    double sigma2;      /* nongenetic (residual) variance */
    double joint_su2;   /* ARG variance component from --joint; 0 if unused.
                           Sequential fitting has no equivalent: there Stage 2
                           estimates it from the RESIDUAL, and that estimate
                           collapses to ~0 (0.008 vs 0.67 joint). */

    /* BayesR-style continuous mixture prior for SNP effects (beta).
     * sigma2_mix[k] is the prior variance for component k; component 0 is
     * always the null (sigma2_mix[0] = 0, no effect).  pi_mix[k] are the
     * mixture proportions, estimated during inference. */
    int     K;           /* number of mixture components */
    double *sigma2_mix;  /* prior variances [K]; owned by sargent_t */
    double *pi_mix;      /* mixture proportions [K]; sum = 1; owned */

    /* Analogous mixture prior for segment weights (ARG random effects).
     * tau2_mix[m] is the prior variance for component m; tau2_mix[0] = 0 (null).
     * phi_mix[m] are the segment mixture proportions, estimated during inference.
     * rho[s*M + m] = q(c_s = m), the per-segment responsibilities (S x M, row-major).
     * ws2[s] = E[w_s^2] = sum_m rho[s*M+m] * tau2_mix[m]. */
    int     M;           /* number of mixture components for segment weights */
    double *tau2_mix;    /* prior variances [M]; tau2_mix[0] = 0; owned */
    double *phi_mix;     /* segment mixture proportions [M]; sum = 1; owned */
    double *rho;         /* segment responsibilities [S * M]; owned */

    /* Sample-rank → individual-index mapping, length 2*n.
     * Built from the tree sequence's individual table at construction time;
     * passed through to tree_cov_matvec so the diploid pairing is always
     * derived from the tree sequence metadata rather than assumed from rank/2. */
    int *hap_to_indiv;

    /* Tree pruning for Stage 2.
     * Adjacent local trees are coarsened into S pruned segments before
     * inference.  S_orig is the original tree count; prune_map[s_orig]
     * gives the pruned segment index for original tree s_orig. */
    int  S_orig;       /* original (unpruned) tree count */
    int *prune_map;    /* prune_map[s_orig] = s_pruned; length S_orig; owned */

    /* Tuning parameters for Stage 2. */
    int nprobe;    /* Hutchinson probe count; default 10 */

    /* Beta(phi_alpha, phi_beta) hyperprior on the activation probability phi
     * for the spike-and-slab post-hoc classification.  Posterior mean update:
     *   phi = (phi_alpha + sum_s rho_s) / (phi_alpha + phi_beta + S)
     * Defaults to 1/1 (uniform).  Set phi_beta >> 1 for a sparse prior. */
    double phi_alpha;
    double phi_beta;

    /* Alternative fitting modes (for systematic comparison).
     * uniform_arg: 1 = single sigma_u^2 shared across all trees in Stage 2
     *   (effectively ARG-BLUP with a uniformly-weighted kinship matrix);
     *   0 = per-tree REML EM (default).
     * em_iters: number of outer iterations of joint variational coordinate
     *   descent in var_joint; 0 = no EM (default; sequential SNP then ARG). */
    int uniform_arg;
    int em_iters;

    /* Block 2.6 — sparse fused-lasso prior on log(ws²_s).
     *   theta_s := log(ws2_s) - mean_log(ws2)
     *   penalty(theta) = lambda1 * sum_s |theta_s|                  (sparsity-on-levels)
     *                  + lambda2 * sum_s |theta_s - theta_{s-1}|    (sparsity-on-differences)
     * Applied as a proximal step at the end of each var_random EM iteration.
     * Defaults (both 0.0) reproduce the unregularised behaviour exactly. */
    double fused_lasso_lambda1;
    double fused_lasso_lambda2;

    /* Precomputed G_all eigendecomposition (optional).
     * When non-NULL and uniform_arg=1, var_random_emma() replaces the
     * CG-based var_random() — reducing Stage-2 from >60 min to <1 s.
     * When non-NULL and uniform_arg=0, used as G_all preconditioner
     * for the per-tree CG solves (cuts iteration count ~64 → ~5).
     * Built by gall_build() or loaded from <prefix>.gall cache;
     * managed externally (sargent_t does not free it). */
    gall_t *gall;
} sargent_t;

sargent_t *sargent_new(geno_reader_t *gr, tsk_treeseq_t *ts, double *y);
void       sargent_free(sargent_t *sg);
