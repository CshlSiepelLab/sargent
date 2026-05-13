#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

#include "sargent.h"
#include "variational.h"
#include "geno.h"

#define VAR_FIXED_MAX_ITER 1000
#define VAR_FIXED_TOL      1e-6

sargent_t *sargent_new(geno_reader_t *gr, tsk_treeseq_t *ts, double *y)
{
    sargent_t *sg = calloc(1, sizeof(*sg));
    if (!sg) return NULL;

    sg->gr = gr;
    sg->ts = ts;
    sg->n  = gr->n_samples;
    sg->p  = gr->n_variants;
    sg->S  = (int)tsk_treeseq_get_num_trees(ts);

    sg->y = malloc((size_t)sg->n * sizeof(double));
    if (!sg->y) { free(sg); return NULL; }
    memcpy(sg->y, y, (size_t)sg->n * sizeof(double));

    /* beta, ws2 remain NULL; scalars zero-initialised by calloc */
    return sg;
}

void sargent_free(sargent_t *sg)
{
    if (!sg) return;
    free(sg->y);
    free(sg->beta);
    free(sg->ws2);
    free(sg);
}

void var_fixed(sargent_t *sg)
{
    int n = sg->n;
    int p = sg->p;

    if (!sg->beta) {
        sg->beta = calloc((size_t)p, sizeof(double));
        if (!sg->beta) {
            fprintf(stderr, "var_fixed: out of memory\n");
            return;
        }
    }

    double *g = malloc((size_t)n * sizeof(double)); /* one marker column */
    double *r = malloc((size_t)n * sizeof(double)); /* running residuals  */
    if (!g || !r) {
        fprintf(stderr, "var_fixed: out of memory\n");
        free(g); free(r);
        return;
    }

    /* Initialise residuals: r = y - X*beta.
       Since beta starts at zero, r = y on the first call. */
    memcpy(r, sg->y, (size_t)n * sizeof(double));
    for (int j = 0; j < p; j++) {
        if (sg->beta[j] == 0.0) continue;
        geno_seek(sg->gr, j);
        geno_read_variant(sg->gr, g, GENO_STANDARDIZE);
        double b = sg->beta[j];
        for (int i = 0; i < n; i++) r[i] -= b * g[i];
    }

    /* Outer loop: sweep until convergence */
    for (int iter = 0; iter < VAR_FIXED_MAX_ITER; iter++) {
        double max_delta = 0.0;

        /* Inner loop: coordinate update over markers */
        geno_rewind(sg->gr);
        for (int j = 0; j < p; j++) {
            if (geno_read_variant(sg->gr, g, GENO_STANDARDIZE) <= 0) break;

            double beta_old = sg->beta[j];

            /* Restore residual for marker j */
            for (int i = 0; i < n; i++) r[i] += beta_old * g[i];

            /* Coordinate update — TODO: replace with variational update */
            double gtr = 0.0;
            for (int i = 0; i < n; i++) gtr += g[i] * r[i];
            double beta_new = gtr / n;

            sg->beta[j] = beta_new;

            /* Remove new contribution from residual */
            for (int i = 0; i < n; i++) r[i] -= beta_new * g[i];

            double d = fabs(beta_new - beta_old);
            if (d > max_delta) max_delta = d;
        }

        /* Update h2_fixed: fraction of phenotypic variance explained */
        double var_fit = 0.0, var_y = 0.0, mean_y = 0.0;
        for (int i = 0; i < n; i++) mean_y += sg->y[i];
        mean_y /= n;
        for (int i = 0; i < n; i++) {
            double fitted = sg->y[i] - r[i];
            var_fit += fitted * fitted;
            double dy = sg->y[i] - mean_y;
            var_y   += dy * dy;
        }
        sg->h2_fixed = (var_y > 0.0) ? var_fit / var_y : 0.0;

        if (max_delta < VAR_FIXED_TOL) {
            fprintf(stderr, "var_fixed: converged at iter %d  h2_fixed=%.4f\n",
                    iter + 1, sg->h2_fixed);
            break;
        }
    }

    free(g);
    free(r);
}
