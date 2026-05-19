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

/* Same as var_fixed but with an explicit additive offset: BayesR is fit as if
   the response were (y - offset).  Used by the joint-EM mode where offset is
   the current random-effect prediction u_hat = K V^{-1} (y - X beta).
   Pass offset == NULL for the default (equivalent to var_fixed). */
void var_fixed_with_offset(sargent_t *sg, const double *offset, int silent);

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

/* Post-convergence pass: compute mixture responsibilities gamma[j*K+k] for
   all p SNPs using the converged sg->beta and sg->sigma2 from var_fixed.
   Returns a malloc'd p*K array (row-major), or NULL on error or if
   sg->beta is NULL.  Caller must free the returned pointer. */
double *var_fixed_posteriors(sargent_t *sg);
