#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

#include "sargent.h"
#include "variational.h"
#include "geno.h"

/* Default BayesR mixture grid: null + three effect-size classes.
 * sigma2_mix values are set relative to var(y) in sargent_new. */
#define VAR_DEFAULT_K       4
#define VAR_FIXED_MAX_ITER  500
#define VAR_FIXED_TOL       1e-6
#define PI_FLOOR            1e-10   /* lower bound on pi_mix to avoid log(0) */

/* -------------------------------------------------------------------------
 * sargent_new / sargent_free
 * ---------------------------------------------------------------------- */

static double vec_var(const double *y, int n)
{
    double mean = 0.0;
    for (int i = 0; i < n; i++) mean += y[i];
    mean /= n;
    double v = 0.0;
    for (int i = 0; i < n; i++) { double d = y[i] - mean; v += d * d; }
    return v / n;
}

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

    /* Initialise residual variance from data */
    double vy = vec_var(y, sg->n);
    sg->sigma2 = vy;

    /* BayesR-style mixture prior: null + three effect-size classes.
     * Scale sigma2_mix by var(y) so the grid is in phenotype units. */
    sg->K = VAR_DEFAULT_K;
    sg->sigma2_mix = malloc((size_t)sg->K * sizeof(double));
    sg->pi_mix     = malloc((size_t)sg->K * sizeof(double));
    if (!sg->sigma2_mix || !sg->pi_mix) {
        sargent_free(sg);
        return NULL;
    }
    sg->sigma2_mix[0] = 0.0;          /* null component */
    sg->sigma2_mix[1] = 1e-4 * vy;
    sg->sigma2_mix[2] = 1e-3 * vy;
    sg->sigma2_mix[3] = 1e-2 * vy;

    /* Start with uniform mixture proportions; updated during inference */
    for (int k = 0; k < sg->K; k++)
        sg->pi_mix[k] = 1.0 / sg->K;

    return sg;
}

void sargent_free(sargent_t *sg)
{
    if (!sg) return;
    free(sg->y);
    free(sg->beta);
    free(sg->ws2);
    free(sg->sigma2_mix);
    free(sg->pi_mix);
    free(sg);
}

/* -------------------------------------------------------------------------
 * var_fixed: Stage 1 variational inference for fixed effects only.
 *
 * Implements the BayesR-style variational coordinate descent with V = σ²I.
 * On each sweep over all p SNPs:
 *
 *   For each j:
 *     u_j = n/σ²                               (constant; standardised x_j)
 *     t_j = x_j^T r_{-j} / σ²
 *     μ_jk = t_j / (u_j + σ_k^{-2})           (posterior mean, component k)
 *     v_jk = 1 / (u_j + σ_k^{-2})             (posterior variance)
 *     γ_jk ∝ π_k N(t_j; 0, u_j σ_k² + 1)     (responsibility)
 *     E[β_j] = Σ_k γ_jk μ_jk
 *
 *   After sweep:
 *     π_k  ← (Σ_j γ_jk) / p
 *     σ²   ← r^T r / n
 * ---------------------------------------------------------------------- */
void var_fixed(sargent_t *sg)
{
    int n = sg->n;
    int p = sg->p;
    int K = sg->K;

    if (!sg->beta) {
        sg->beta = calloc((size_t)p, sizeof(double));
        if (!sg->beta) { fprintf(stderr, "var_fixed: out of memory\n"); return; }
    }

    double *g         = malloc((size_t)n * sizeof(double));
    double *r         = malloc((size_t)n * sizeof(double));
    double *log_g     = malloc((size_t)K * sizeof(double));
    double *gamma     = malloc((size_t)K * sizeof(double));
    double *gamma_sum = malloc((size_t)K * sizeof(double));
    if (!g || !r || !log_g || !gamma || !gamma_sum) {
        fprintf(stderr, "var_fixed: out of memory\n");
        free(g); free(r); free(log_g); free(gamma); free(gamma_sum);
        return;
    }

    /* Initialise residuals r = y - X β.
     * On first call β = 0 so r = y; on re-entry the current β is used. */
    memcpy(r, sg->y, (size_t)n * sizeof(double));
    for (int j = 0; j < p; j++) {
        if (sg->beta[j] == 0.0) continue;
        geno_seek(sg->gr, j);
        geno_read_variant(sg->gr, g, GENO_STANDARDIZE);
        double b = sg->beta[j];
        for (int i = 0; i < n; i++) r[i] -= b * g[i];
    }

    double var_y = vec_var(sg->y, n);

    for (int iter = 0; iter < VAR_FIXED_MAX_ITER; iter++) {

        /* u_bar = n/σ² is constant across all SNPs for standardised genotypes */
        double u_bar = (double)n / sg->sigma2;

        double max_delta = 0.0;
        memset(gamma_sum, 0, (size_t)K * sizeof(double));
        int n_read = 0;

        geno_rewind(sg->gr);
        for (int j = 0; j < p; j++) {
            if (geno_read_variant(sg->gr, g, GENO_STANDARDIZE) <= 0) break;
            n_read++;

            double beta_old = sg->beta[j];

            /* Restore leave-one-out residual: r_{-j} = r + β_j x_j */
            for (int i = 0; i < n; i++) r[i] += beta_old * g[i];

            /* t_j = x_j^T V^{-1} r_{-j} = (1/σ²) x_j^T r_{-j} */
            double xr = 0.0;
            for (int i = 0; i < n; i++) xr += g[i] * r[i];
            double t_j = xr / sg->sigma2;

            /* Log unnormalised responsibilities.
             * For component k: log γ_jk ∝ log π_k + log N(t_j; 0, u_bar σ_k² + 1).
             * Constant -0.5 log(2π) is dropped (cancels in normalisation).
             * Null component (σ_k² = 0): log N(t_j; 0, 1) = -0.5 t_j². */
            double lmax = -1e300;
            for (int k = 0; k < K; k++) {
                double sk2 = sg->sigma2_mix[k];
                double log_pi_k = log(sg->pi_mix[k] + PI_FLOOR);
                if (sk2 == 0.0) {
                    log_g[k] = log_pi_k - 0.5 * t_j * t_j;
                } else {
                    double s_k = u_bar * sk2 + 1.0;
                    log_g[k] = log_pi_k - 0.5 * (log(s_k) + t_j * t_j / s_k);
                }
                if (log_g[k] > lmax) lmax = log_g[k];
            }

            /* Softmax normalisation (log-sum-exp) */
            double sum_g = 0.0;
            for (int k = 0; k < K; k++) {
                gamma[k] = exp(log_g[k] - lmax);
                sum_g += gamma[k];
            }
            for (int k = 0; k < K; k++) gamma[k] /= sum_g;

            /* Posterior mean: E[β_j] = Σ_k γ_jk μ_jk
             * μ_jk = t_j v_jk = t_j / (u_bar + 1/σ_k²);  μ_j0 = 0 (null). */
            double beta_new = 0.0;
            for (int k = 0; k < K; k++) {
                if (gamma[k] == 0.0 || sg->sigma2_mix[k] == 0.0) continue;
                double v_jk   = 1.0 / (u_bar + 1.0 / sg->sigma2_mix[k]);
                double mu_jk  = t_j * v_jk;
                beta_new     += gamma[k] * mu_jk;
            }

            sg->beta[j] = beta_new;

            /* Update residual: subtract new contribution */
            for (int i = 0; i < n; i++) r[i] -= beta_new * g[i];

            /* Accumulate responsibilities for π update */
            for (int k = 0; k < K; k++) gamma_sum[k] += gamma[k];

            double delta = fabs(beta_new - beta_old);
            if (delta > max_delta) max_delta = delta;
        }

        /* Update mixture proportions: π_k = (Σ_j γ_jk) / p */
        if (n_read > 0) {
            for (int k = 0; k < K; k++)
                sg->pi_mix[k] = gamma_sum[k] / n_read;
        }

        /* Update residual variance: σ² = r^T r / n */
        double rr = 0.0;
        for (int i = 0; i < n; i++) rr += r[i] * r[i];
        sg->sigma2 = rr / n;
        if (sg->sigma2 < 1e-10 * var_y)
            sg->sigma2 = 1e-10 * var_y;  /* floor to avoid u_bar → ∞ */

        /* Update h2_fixed: fraction of phenotypic variance explained by X β̂ */
        sg->h2_fixed  = (var_y > 0.0) ? 1.0 - sg->sigma2 / var_y : 0.0;

        fprintf(stderr, "var_fixed iter %3d: max_delta=%.2e  sigma2=%.4f"
                "  h2_fixed=%.4f  pi=[",
                iter + 1, max_delta, sg->sigma2, sg->h2_fixed);
        for (int k = 0; k < K; k++)
            fprintf(stderr, " %.3f", sg->pi_mix[k]);
        fprintf(stderr, " ]\n");

        if (max_delta < VAR_FIXED_TOL) break;
    }

    free(g);
    free(r);
    free(log_g);
    free(gamma);
    free(gamma_sum);
}
