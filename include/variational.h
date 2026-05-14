#pragma once

#include "sargent.h"

/* Stage 1: infer fixed effects (beta) via BayesR-style variational coordinate
   descent with V = sigma^2 I.  Allocates sg->beta if not already set.
   If silent != 0, per-iteration progress is suppressed. */
void var_fixed(sargent_t *sg, int silent);

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
