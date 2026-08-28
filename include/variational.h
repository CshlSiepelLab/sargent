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

#include "sargent.h"
#include "gall.h"

/* Stage 1: infer fixed effects (beta) via BayesR-style variational coordinate
   descent with V = sigma^2 I.  Allocates sg->beta if not already set.
   If silent != 0, per-iteration progress is suppressed. */
void var_fixed(sargent_t *sg, int silent);

/* Select the Stage-1 inference engine.  on=0 (default) keeps the variational
 * fit; on=1 replaces it with a BayesR Gibbs sampler (SSVS), reporting the
 * post-burn-in posterior mean of beta.  niter<=0 / burn<0 keep the defaults
 * (3000 / 1000).  Set from --stage1-gibbs, --gibbs-iter, --gibbs-burn. */
void var_fixed_set_gibbs(int on, int niter, int burn);

/* JOINT estimation of beta and the ARG variance components in one Gibbs chain,
 * in the eigenbasis cached in .gall.  Replaces the sequential scheme, which
 * gives the ARG nothing: sargent_u equalled bayesR to six decimals at every phi
 * level while arg_lmm alone beat bayesR by up to 0.05.  Measured +0.0218 r2
 * over sequential (19/21 cells) at 0.83x the cost.  Needs a .gall whose n
 * matches; returns 0 on success. */
int var_joint_gibbs(sargent_t *sg, const gall_t *gl, int silent,
                    int niter, int burn, uint64_t seed,
                    double *su2_out, double *sigma2_out);

/* Joint variational coordinate descent: each outer iteration performs one
   V-solve, one sweep of beta updates using the V-aware precision tr(V^{-1}),
   and one REML EM step on (ws^2, sigma^2).  Initialised with beta=0 and a
   warm ws^2 so the first beta sweep already sees a non-trivial random offset.
   Iterates sg->em_iters times; post-hoc rho computed at end. */
void var_joint(sargent_t *sg, int silent);

/* Stage 2: infer segment weights (ws2) and sigma^2 via REML-style score
   updates, using CG for whitened residuals and Hutchinson trace estimates.
   Fixed effects are held at their Stage 1 values.
   Allocates sg->ws2, sg->rho if not already set.
   If silent != 0, per-iteration progress is suppressed. */
void var_random(sargent_t *sg, int silent);

/* Stage 2 for uniform-ARG mode using EMMA profile-likelihood search.
 * Requires sg->gall to be set (non-NULL) and sg->uniform_arg = 1.
 * Replaces var_random() in that configuration: no CG solves, no Hutchinson
 * traces — just one O(n²) transform plus a 100-step golden-section search.
 * On return, sg->ws2[s] = ws2_total/S for all s, and sg->sigma2, sg->h2_arg
 * are set consistently.
 * If sg->gall is NULL, falls back to var_random(). */
void var_random_emma(sargent_t *sg, int silent);

/* Post-convergence pass: compute mixture responsibilities gamma[j*K+k] for
   all p SNPs using the converged sg->beta and sg->sigma2 from var_fixed.
   Returns a malloc'd p*K array (row-major), or NULL on error or if
   sg->beta is NULL.  Caller must free the returned pointer. */
double *var_fixed_posteriors(sargent_t *sg);
