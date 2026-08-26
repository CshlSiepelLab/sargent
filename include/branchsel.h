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

#include <stdint.h>
#include <tskit.h>

/* Branch selection — the Stage-2 replacement.
 *
 * Ranks every ARG branch by the marginal LMM score statistic against the
 * uniform-ARG background,
 *
 *     z_b = | x_b' P y | / sqrt( x_b' P x_b )
 *
 * where x_b is the per-individual carrier count of branch b (0/1/2 for a
 * diploid) and P is the REML projector of the fitted uniform-ARG model,
 *
 *     P = V^-1 - V^-1 1 1' V^-1 / (1' V^-1 1),   V = sigma_u^2 K_uni + sigma_e^2 I.
 *
 * The top k branches then form a second kinship K_sel = Z_sel Z_sel' / k, and a
 * two-component REML (K_uni + K_sel) replaces the per-tree EM.  z_b is invariant
 * to centering and scaling of x_b (P 1 = 0 kills the mean and the scale cancels
 * between numerator and denominator), so the raw carrier counts are used here
 * and standardization is applied only when building Z_sel.
 *
 * Both the score and the selected kinship are computed on ONE tree sequence —
 * the full (train + test) ARG — with the score restricted to train rows.  This
 * is deliberate: a branch selected in the train-only ARG has no well-defined
 * carrier set among test individuals, because simplification collapses a whole
 * root-ward path of full-ARG branches (which differ in their test carriers)
 * into a single train-ARG branch.  The test individuals' phenotypes are never
 * touched; only their position in the genealogy is used, exactly as for the
 * uniform kinship.
 */

typedef struct {
    int k;          /* branches to select; <= 0 -> BRANCHSEL_DEFAULT_K (25) */
    int max_trees;  /* scan this many evenly-spaced local trees; 0 = every tree.
                       The validated Python prototype used 200. */
    int silent;
} branchsel_opts_t;

#define BRANCHSEL_DEFAULT_K      25
#define BRANCHSEL_DEFAULT_TREES 200

typedef struct {
    int     n_total;    /* individuals in the tree sequence (train + test) */
    int     n_train;    /* leading rows used for scoring */
    int     n_branch;   /* distinct, non-degenerate branches scored */
    int     n_trees;    /* local trees actually scanned */
    int     k;          /* branches selected (min(opts.k, n_branch)) */
    double *Z;          /* [n_total * k] row-major; column t is branch t
                           standardized by its TRAIN mean and population sd */
    double *score;      /* [k] the |z| of each selected branch, descending */
    double  seconds;
} branchsel_t;

/* Score every branch and select the top k.
 *
 * ts:           full ARG; individuals [0, n_train) are train, the rest test.
 * hap_to_indiv: from tsk_build_hap_to_indiv(ts).
 * P:            [n_train * n_train] REML projector, row-major symmetric.
 * y_centered:   [n_train] train phenotype minus its mean.
 *
 * Returns 0 on success (out is filled and owned by the caller, free with
 * branchsel_free), -1 on error. */
int branchsel_run(const tsk_treeseq_t *ts,
                  const int           *hap_to_indiv,
                  int                  n_total,
                  int                  n_train,
                  const double        *P,
                  const double        *y_centered,
                  const branchsel_opts_t *opts,
                  branchsel_t         *out);

void branchsel_free(branchsel_t *r);

/* K_sel blocks from Z: K = Z Z' / k.
 * K_train is [n_train * n_train], K_cross is [(n_total-n_train) * n_train],
 * both row-major, matching predict.c's kinship-block convention. */
void branchsel_kinship(const branchsel_t *r,
                       double *K_train,
                       double *K_cross);

/* FNV-1a 64 over the n_total carrier bytes of a branch — the fingerprint written
 * to the SARGENT_BRANCHSEL_DUMP file so an external reference implementation can
 * align its branches with ours.  Exposed for the unit test. */
uint64_t branchsel_fingerprint(const unsigned char *carriers, int n);
