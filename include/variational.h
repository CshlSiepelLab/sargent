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

/* Stage 1: infer fixed effects (beta) via BayesR-style variational coordinate
   descent with V = sigma^2 I.  Allocates sg->beta if not already set.
   If silent != 0, per-iteration progress is suppressed. */
void var_fixed(sargent_t *sg, int silent);

/* Joint variational coordinate descent: each outer iteration performs one
   V-solve, one sweep of beta updates using the V-aware precision tr(V^{-1}),
   and one REML EM step on (ws^2, sigma^2).  Initialised with beta=0 and a
   warm ws^2 so the first beta sweep already sees a non-trivial random offset.
   Iterates sg->em_iters times; post-hoc rho computed at end. */
void var_joint(sargent_t *sg, int silent);

/* Method of moments variance-component estimation.
 *
 * Fits the three-VC model
 *    y = X beta + u + epsilon,   beta ~ N(0, sigma2_beta_tilde / p · I_p),
 *    u ~ N(0, sigma2_g K_g),     epsilon ~ N(0, sigma2 I),
 * where K_X = X X^T / p and K_g = (1/S) sum_s G_s.
 *
 * Solves the 3-equation linear system from E[y^T M y] = tr(M V) for
 * M in {I, K_X, K_g} via Hutchinson stochastic trace estimates.
 * Negative variance estimates are clipped to zero.
 *
 * Stores results in sg->h2_fixed, sg->h2_random, sg->h2, sg->sigma2.
 * Unlike var_fixed + var_random or var_joint, the MoM estimator does not
 * suffer from the fixed/random identifiability bias: it estimates
 * variance components without first identifying beta or per-tree ws2,
 * so cannibalisation does not affect it.
 *
 * Does not populate sg->beta or sg->ws2 — those remain whatever the
 * caller initialised them to (typically NULL on first invocation). */
void var_mom(sargent_t *sg, int silent);

/* Stage 2: infer segment weights (ws2) and sigma^2 via REML-style score
   updates, using CG for whitened residuals and Hutchinson trace estimates.
   Fixed effects are held at their Stage 1 values.
   Allocates sg->ws2, sg->rho if not already set.
   If silent != 0, per-iteration progress is suppressed. */
void var_random(sargent_t *sg, int silent);

/* Post-convergence pass: compute mixture responsibilities gamma[j*K+k] for
   all p SNPs using the converged sg->beta and sg->sigma2 from var_fixed.
   Returns a malloc'd p*K array (row-major), or NULL on error or if
   sg->beta is NULL.  Caller must free the returned pointer. */
double *var_fixed_posteriors(sargent_t *sg);
