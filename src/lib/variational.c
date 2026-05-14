/*
 * SARGENT: Sparse ARG-based ENcoding of Traits
 *
 * Copyright (c) 2026, Cold Spring Harbor Laboratory
 * All rights reserved.
 *
 * This file is part of SARGENT and is distributed under the BSD 3-Clause License.
 * See the LICENSE file in the project root for details.
 */

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

#include "sargent.h"
#include "variational.h"
#include "geno.h"
#include "conjgrad.h"
#include "hutchinson.h"

/* Default BayesR mixture grid: null + three effect-size classes.
 * sigma2_mix values are set relative to var(y) in sargent_new. */
#define VAR_DEFAULT_K        4
#define VAR_DEFAULT_M        4
#define VAR_FIXED_MAX_ITER   500
#define VAR_RANDOM_MAX_ITER  200
#define VAR_FIXED_TOL        1e-6
#define VAR_RANDOM_TOL       1e-4
#define PI_FLOOR             1e-10   /* lower bound on pi/phi_mix to avoid log(0) */
#define SCORE_CAP            50.0    /* clip tau^2*(A-B)/2 before softmax */

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

    /* Segment-weight mixture prior.
     * tau2_mix is scaled by vy/S so that the total random-effect variance
     * (if all S segments are in the largest component) is of order vy. */
    sg->M = VAR_DEFAULT_M;
    sg->tau2_mix = malloc((size_t)sg->M * sizeof(double));
    sg->phi_mix  = malloc((size_t)sg->M * sizeof(double));
    sg->rho      = malloc((size_t)sg->S * sg->M * sizeof(double));
    if (!sg->tau2_mix || !sg->phi_mix || !sg->rho) {
        sargent_free(sg);
        return NULL;
    }
    double seg_scale = vy / (sg->S > 0 ? sg->S : 1);
    sg->tau2_mix[0] = 0.0;
    sg->tau2_mix[1] = 0.01  * seg_scale;
    sg->tau2_mix[2] = 0.1   * seg_scale;
    sg->tau2_mix[3] = 1.0   * seg_scale;

    /* Uniform initial responsibilities and proportions */
    for (int m = 0; m < sg->M; m++)
        sg->phi_mix[m] = 1.0 / sg->M;
    for (int s = 0; s < sg->S; s++)
        for (int m = 0; m < sg->M; m++)
            sg->rho[s * sg->M + m] = 1.0 / sg->M;

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
    free(sg->tau2_mix);
    free(sg->phi_mix);
    free(sg->rho);
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
void var_fixed(sargent_t *sg, int silent)
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

        if (!silent) {
            fprintf(stderr, "var_fixed iter %3d: max_delta=%.2e  sigma2=%.4f"
                    "  h2_fixed=%.4f  pi=[",
                    iter + 1, max_delta, sg->sigma2, sg->h2_fixed);
            for (int k = 0; k < K; k++)
                fprintf(stderr, " %.3f", sg->pi_mix[k]);
            fprintf(stderr, " ]\n");
        }

        if (max_delta < VAR_FIXED_TOL) break;
    }

    free(g);
    free(r);
    free(log_g);
    free(gamma);
    free(gamma_sum);
}

/* =========================================================================
 * Stage 2: var_random — segment weight inference
 * ======================================================================= */

/* V matvec: out = (Σ_s ws2[s] G_s + sigma2 I) * in.
 * Makes one sequential pass over the tree sequence. */
typedef struct {
    const tsk_treeseq_t *ts;
    const double        *ws2;    /* E[w_s^2] for each tree, length S */
    double               sigma2;
    int                  n;
    int                  S;
    double              *work;      /* 2 * n_nodes scratch (borrowed) */
    tsk_id_t            *order_buf; /* n_nodes scratch (borrowed) */
    double              *gv;        /* n scratch (borrowed) */
} v_matvec_data_t;

static void v_matvec(double *out, const double *in, int n, void *data)
{
    v_matvec_data_t *d = data;

    for (int i = 0; i < n; i++) out[i] = d->sigma2 * in[i];

    tsk_tree_t tree;
    int ret = tsk_tree_init(&tree, d->ts, TSK_NO_SAMPLE_COUNTS);
    if (ret < 0) {
        fprintf(stderr, "v_matvec: tsk_tree_init: %s\n", tsk_strerror(ret));
        return;
    }
    int s = 0;
    for (ret = tsk_tree_first(&tree); ret == TSK_TREE_OK && s < d->S;
         ret = tsk_tree_next(&tree), s++) {
        if (d->ws2[s] == 0.0) continue;
        tree_cov_matvec(d->gv, in, d->n, &tree, d->work, d->order_buf);
        for (int i = 0; i < n; i++) out[i] += d->ws2[s] * d->gv[i];
    }
    if (ret < 0)
        fprintf(stderr, "v_matvec: tree error: %s\n", tsk_strerror(ret));
    tsk_tree_free(&tree);
}

/* Compute A[s] = r_tilde^T G_s r_tilde for all S trees in one pass. */
static void compute_quad_all(const tsk_treeseq_t *ts, int n, int S,
                              const double *r_tilde,
                              double *A,
                              double *work, tsk_id_t *order_buf, double *gv)
{
    tsk_tree_t tree;
    int ret = tsk_tree_init(&tree, ts, TSK_NO_SAMPLE_COUNTS);
    if (ret < 0) {
        fprintf(stderr, "compute_quad_all: tsk_tree_init: %s\n", tsk_strerror(ret));
        for (int s = 0; s < S; s++) A[s] = 0.0;
        return;
    }
    int s = 0;
    for (ret = tsk_tree_first(&tree); ret == TSK_TREE_OK && s < S;
         ret = tsk_tree_next(&tree), s++) {
        tree_cov_matvec(gv, r_tilde, n, &tree, work, order_buf);
        double a = 0.0;
        for (int i = 0; i < n; i++) a += r_tilde[i] * gv[i];
        A[s] = a;
    }
    /* zero any remaining entries if tree sequence shorter than expected */
    for (; s < S; s++) A[s] = 0.0;
    if (ret < 0)
        fprintf(stderr, "compute_quad_all: tree error: %s\n", tsk_strerror(ret));
    tsk_tree_free(&tree);
}

void var_random(sargent_t *sg, int silent)
{
    int n = sg->n;
    int S = sg->S;
    int M = sg->M;

    /* Allocate ws2 if needed; start from zero (will be overwritten on iter 0) */
    if (!sg->ws2) {
        sg->ws2 = calloc((size_t)S, sizeof(double));
        if (!sg->ws2) { fprintf(stderr, "var_random: out of memory\n"); return; }
    }

    /* Working arrays */
    double   *r        = malloc((size_t)n * sizeof(double));
    double   *r_tilde  = calloc((size_t)n, sizeof(double)); /* warm-start: keep between iters */
    double   *A        = malloc((size_t)S * sizeof(double));
    double   *B        = malloc((size_t)S * sizeof(double));
    double   *log_rho  = malloc((size_t)M * sizeof(double));
    double   *phi_sum  = malloc((size_t)M * sizeof(double));

    tsk_size_t n_nodes  = tsk_treeseq_get_num_nodes(sg->ts);
    double    *work     = malloc(2 * n_nodes * sizeof(double));
    tsk_id_t  *order_buf = malloc(n_nodes * sizeof(tsk_id_t));
    double    *gv        = malloc((size_t)n * sizeof(double));

    cg_state_t *cg = cg_state_new(n, NULL);

    if (!r || !r_tilde || !A || !B || !log_rho || !phi_sum ||
        !work || !order_buf || !gv || !cg) {
        fprintf(stderr, "var_random: out of memory\n");
        goto cleanup;
    }

    /* Build residuals r = y - X E[beta] (fixed effects held constant) */
    memcpy(r, sg->y, (size_t)n * sizeof(double));
    if (sg->beta) {
        double *g = malloc((size_t)n * sizeof(double));
        if (!g) { fprintf(stderr, "var_random: out of memory\n"); goto cleanup; }
        geno_rewind(sg->gr);
        for (int j = 0; j < sg->p; j++) {
            if (geno_read_variant(sg->gr, g, GENO_STANDARDIZE) <= 0) break;
            if (sg->beta[j] == 0.0) continue;
            for (int i = 0; i < n; i++) r[i] -= sg->beta[j] * g[i];
        }
        free(g);
    }

    double var_y = vec_var(sg->y, n);

    /* V matvec callback data — ws2 pointer is live (updated in place) */
    v_matvec_data_t vd = {
        .ts = sg->ts, .ws2 = sg->ws2, .sigma2 = sg->sigma2,
        .n = n, .S = S,
        .work = work, .order_buf = order_buf, .gv = gv,
    };

    hutch_opts_t hopts = {
        .nprobe = 50, .batch_size = 0, .use_gaussian = 0, .seed = 0,
    };

    double old_sigma2 = sg->sigma2;

    for (int iter = 0; iter < VAR_RANDOM_MAX_ITER; iter++) {

        /* Keep V matvec in sync with current sigma2 (ws2 pointer stays live) */
        vd.sigma2 = sg->sigma2;

        /* --- a. Solve V r_tilde = r (warm-start from previous r_tilde) --- */
        cg_solve(cg, v_matvec, &vd, r, r_tilde);

        /* --- b. A_0 = r_tilde^T r_tilde ----------------------------------- */
        double A_0 = 0.0;
        for (int i = 0; i < n; i++) A_0 += r_tilde[i] * r_tilde[i];

        /* --- c. A_s = r_tilde^T G_s r_tilde for all s -------------------- */
        compute_quad_all(sg->ts, n, S, r_tilde, A, work, order_buf, gv);

        /* --- d. B_s = tr(V^{-1} G_s) and B_0 = tr(V^{-1}) via Hutchinson - */
        double B_0 = 0.0;
        hutch_trace_all(cg, v_matvec, &vd, sg->ts, n, S, B, &B_0, &hopts);

        /* --- e. Update rho_sm, ws2[s], phi_mix ---------------------------- */
        memset(phi_sum, 0, (size_t)M * sizeof(double));
        double max_delta_ws2 = 0.0;

        for (int s = 0; s < S; s++) {
            double lmax = -1e300;
            for (int m = 0; m < M; m++) {
                double score = sg->tau2_mix[m] * 0.5 * (A[s] - B[s]);
                /* Clip to prevent exp overflow; large positive score drives
                 * all weight to component m, which is the correct behavior. */
                score = fmax(-SCORE_CAP, fmin(SCORE_CAP, score));
                log_rho[m] = log(sg->phi_mix[m] + PI_FLOOR) + score;
                if (log_rho[m] > lmax) lmax = log_rho[m];
            }

            /* Softmax normalisation */
            double sum = 0.0;
            for (int m = 0; m < M; m++) {
                sg->rho[s * M + m] = exp(log_rho[m] - lmax);
                sum += sg->rho[s * M + m];
            }
            double ws2_old = sg->ws2[s];
            double ws2_new = 0.0;
            for (int m = 0; m < M; m++) {
                sg->rho[s * M + m] /= sum;
                ws2_new += sg->rho[s * M + m] * sg->tau2_mix[m];
                phi_sum[m] += sg->rho[s * M + m];
            }
            sg->ws2[s] = ws2_new;
            double d = fabs(ws2_new - ws2_old);
            if (d > max_delta_ws2) max_delta_ws2 = d;
        }

        /* Update segment mixture proportions: phi_m = (Σ_s rho_sm) / S */
        for (int m = 0; m < M; m++)
            sg->phi_mix[m] = phi_sum[m] / S;

        /* --- f. Update sigma2: sigma2 <- sigma2 * A_0 / B_0 (REML step) -- */
        if (B_0 > 0.0) {
            sg->sigma2 *= A_0 / B_0;
            if (sg->sigma2 < 1e-10 * var_y) sg->sigma2 = 1e-10 * var_y;
            if (sg->sigma2 > var_y)         sg->sigma2 = var_y;
        }

        /* Update summary statistics */
        sg->h2_random = (var_y > 0.0) ? fmax(0.0, 1.0 - sg->sigma2/var_y - sg->h2_fixed) : 0.0;
        sg->h2        = sg->h2_fixed + sg->h2_random;

        double rel_change = fabs(sg->sigma2 - old_sigma2) / (old_sigma2 + 1e-15);

        if (!silent) {
            fprintf(stderr, "var_random iter %3d: sigma2=%.4f  h2_rand=%.4f"
                    "  max_dws2=%.2e  B0=%.2e  rel_dsig=%.2e  phi=[",
                    iter + 1, sg->sigma2, sg->h2_random,
                    max_delta_ws2, B_0, rel_change);
            for (int m = 0; m < M; m++)
                fprintf(stderr, " %.3f", sg->phi_mix[m]);
            fprintf(stderr, " ]\n");
        }

        old_sigma2 = sg->sigma2;
        if (rel_change < VAR_RANDOM_TOL && max_delta_ws2 < VAR_RANDOM_TOL)
            break;
    }

cleanup:
    cg_state_free(cg);
    free(r);
    free(r_tilde);
    free(A);
    free(B);
    free(log_rho);
    free(phi_sum);
    free(work);
    free(order_buf);
    free(gv);
}

/* =========================================================================
 * var_fixed_posteriors: post-convergence gamma_jk computation
 * ======================================================================= */

double *var_fixed_posteriors(sargent_t *sg)
{
    if (!sg->beta) return NULL;

    int n = sg->n;
    int p = sg->p;
    int K = sg->K;

    double *gamma_all = malloc((size_t)p * K * sizeof(double));
    double *g         = malloc((size_t)n * sizeof(double));
    double *r         = malloc((size_t)n * sizeof(double));
    double *log_g     = malloc((size_t)K * sizeof(double));
    double *gamma     = malloc((size_t)K * sizeof(double));
    if (!gamma_all || !g || !r || !log_g || !gamma) {
        fprintf(stderr, "var_fixed_posteriors: out of memory\n");
        free(gamma_all); free(g); free(r); free(log_g); free(gamma);
        return NULL;
    }

    /* Rebuild full residuals r = y - X beta */
    memcpy(r, sg->y, (size_t)n * sizeof(double));
    geno_rewind(sg->gr);
    for (int j = 0; j < p; j++) {
        if (geno_read_variant(sg->gr, g, GENO_STANDARDIZE) <= 0) break;
        if (sg->beta[j] == 0.0) continue;
        for (int i = 0; i < n; i++) r[i] -= sg->beta[j] * g[i];
    }

    double u_bar = (double)n / sg->sigma2;

    /* One sequential pass: for each SNP compute leave-one-out t_j and gamma_jk */
    geno_rewind(sg->gr);
    for (int j = 0; j < p; j++) {
        if (geno_read_variant(sg->gr, g, GENO_STANDARDIZE) <= 0) break;

        /* r_{-j} = r + beta_j x_j */
        for (int i = 0; i < n; i++) r[i] += sg->beta[j] * g[i];

        double xr = 0.0;
        for (int i = 0; i < n; i++) xr += g[i] * r[i];
        double t_j = xr / sg->sigma2;

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
        double sum_g = 0.0;
        for (int k = 0; k < K; k++) {
            gamma[k] = exp(log_g[k] - lmax);
            sum_g += gamma[k];
        }
        for (int k = 0; k < K; k++)
            gamma_all[j * K + k] = gamma[k] / sum_g;

        /* Restore: r -= beta_j x_j */
        for (int i = 0; i < n; i++) r[i] -= sg->beta[j] * g[i];
    }

    free(g); free(r); free(log_g); free(gamma);
    return gamma_all;
}
