#pragma once

#include <tskit.h>
#include "geno.h"

typedef struct {
    int n;              /* number of individuals */
    int p;              /* number of markers */
    int S;              /* number of segments (local trees) */
    geno_reader_t  *gr; /* genotype reader (borrowed) */
    tsk_treeseq_t  *ts; /* tree sequence (borrowed) */
    double         *y;  /* phenotype vector [n] */
    double        *beta;   /* fixed effects [p], NULL until inferred */
    double        *ws2;    /* segment weights [S], NULL until inferred */
    double h2;          /* total heritability */
    double h2_fixed;    /* fixed-effect heritability component */
    double h2_random;   /* random-effect heritability component */
    double sigma2;      /* nongenetic (residual) variance */

    /* BayesR-style continuous mixture prior for fixed effects (beta).
     * sigma2_mix[k] is the prior variance for component k; component 0 is
     * always the null (sigma2_mix[0] = 0, no effect).  pi_mix[k] are the
     * mixture proportions, estimated during inference. */
    int     K;           /* number of mixture components */
    double *sigma2_mix;  /* prior variances [K]; owned by sargent_t */
    double *pi_mix;      /* mixture proportions [K]; sum = 1; owned */

    /* Analogous mixture prior for segment weights (random effects).
     * tau2_mix[m] is the prior variance for component m; tau2_mix[0] = 0 (null).
     * phi_mix[m] are the segment mixture proportions, estimated during inference.
     * rho[s*M + m] = q(c_s = m), the per-segment responsibilities (S x M, row-major).
     * ws2[s] = E[w_s^2] = sum_m rho[s*M+m] * tau2_mix[m]. */
    int     M;           /* number of mixture components for segment weights */
    double *tau2_mix;    /* prior variances [M]; tau2_mix[0] = 0; owned */
    double *phi_mix;     /* segment mixture proportions [M]; sum = 1; owned */
    double *rho;         /* segment responsibilities [S * M]; owned */
} sargent_t;

sargent_t *sargent_new(geno_reader_t *gr, tsk_treeseq_t *ts, double *y);
void       sargent_free(sargent_t *sg);
