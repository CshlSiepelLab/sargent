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
#define VAR_RANDOM_MAX_ITER  500
#define VAR_FIXED_TOL        1e-4
#define VAR_RANDOM_TOL       1e-4   /* convergence threshold for rel_dsig */
#define PI_FLOOR             1e-10   /* lower bound on pi/phi_mix to avoid log(0) */
#define SCORE_CAP            50.0    /* clip tau^2*(A-B)/2 before softmax */
#define VAR_N_SEGMENTS       50      /* target pruned segment count for Stage 2 */

/* Map each of the S_orig original trees to one of n_target pruned segments
 * using equal-count binning: segment sp = floor(s * n_target / S_orig). */
static int *build_prune_map(int S_orig, int n_target)
{
    if (S_orig <= 0) {
        int *pm = malloc(sizeof(int));
        return pm;
    }
    int *pm = malloc((size_t)S_orig * sizeof(int));
    if (!pm) return NULL;
    if (n_target <= 0 || S_orig <= n_target) {
        for (int s = 0; s < S_orig; s++) pm[s] = s;
    } else {
        for (int s = 0; s < S_orig; s++)
            pm[s] = (int)((long long)s * n_target / S_orig);
    }
    return pm;
}

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
    sg->S_orig = (int)tsk_treeseq_get_num_trees(ts);

    sg->prune_map = build_prune_map(sg->S_orig, VAR_N_SEGMENTS);
    if (!sg->prune_map) { free(sg); return NULL; }
    int S_pruned = (sg->S_orig > 0) ? sg->prune_map[sg->S_orig - 1] + 1 : 0;
    sg->S = S_pruned;

    if (sg->S_orig != sg->S)
        fprintf(stderr, "sargent: pruned %d trees → %d segments\n",
                sg->S_orig, sg->S);

    sg->y = malloc((size_t)sg->n * sizeof(double));
    if (!sg->y) { free(sg); return NULL; }
    memcpy(sg->y, y, (size_t)sg->n * sizeof(double));

    /* Center phenotype: LMM/BayesR are translation-invariant only when y is
     * zero-mean.  Subtract sample mean so the intercept is implicit. */
    double ymean = 0.0;
    for (int i = 0; i < sg->n; i++) ymean += sg->y[i];
    ymean /= sg->n;
    for (int i = 0; i < sg->n; i++) sg->y[i] -= ymean;

    /* Initialise residual variance from data */
    double vy = vec_var(sg->y, sg->n);
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
     * tau2_mix is scaled so that the total random-effect variance
     * (if all S segments are in the largest component) is of order vy.
     *
     * G_s entries are in branch-length units (e.g. generations), so we
     * calibrate by estimating E[tr(G_s)/n] via a single Hutchinson probe.
     * This makes tau2_mix independent of the branch-length scale. */
    /* Build diploid mapping first — needed for both G_s tree walks below. */
    sg->hap_to_indiv = tsk_build_hap_to_indiv(ts);
    if (!sg->hap_to_indiv) {
        sargent_free(sg);
        return NULL;
    }

    sg->M = VAR_DEFAULT_M;
    sg->tau2_mix = malloc((size_t)sg->M * sizeof(double));
    sg->phi_mix  = malloc((size_t)sg->M * sizeof(double));
    sg->rho      = malloc((size_t)sg->S * sg->M * sizeof(double));
    if (!sg->tau2_mix || !sg->phi_mix || !sg->rho) {
        sargent_free(sg);
        return NULL;
    }
    /* Estimate E[tr(G_s)/n] so tau2_mix is in phenotype-variance units
     * regardless of branch-length scale (e.g. generations vs coalescent). */
    double gscale    = estimate_g_scale(ts, sg->hap_to_indiv, sg->n, 0u);
    /* Calibrate tau2_mix so total initial random effect ≈ vy.
     * Use S_orig (not S_pruned) because v_matvec sums ALL original trees,
     * each contributing ws2[prune_map[s]] * G_s. */
    double seg_scale = vy / ((double)(sg->S_orig > 0 ? sg->S_orig : 1) * fmax(gscale, 1.0));
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

    sg->nprobe = 10;   /* default; override via sg->nprobe before calling var_random */

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
    free(sg->hap_to_indiv);
    free(sg->prune_map);
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
    double old_sigma2 = sg->sigma2;

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

            /* Log unnormalised responsibilities (Wakefield ABF form).
             * z_j = t_j / sqrt(u_bar) is the standardised z-score (~N(0,1) under null).
             * Null  (sk²=0): log γ ∝ log π_k − ½ z_j²
             * Other (sk²>0): log γ ∝ log π_k − ½(log s_k + z_j²/s_k)
             * where s_k = u_bar·σ_k² + 1.  Constant -½ log(2π) cancels. */
            double z_sq = t_j * t_j / u_bar;  /* ~chi2(1) under null */
            double lmax = -1e300;
            for (int k = 0; k < K; k++) {
                double sk2 = sg->sigma2_mix[k];
                double log_pi_k = log(sg->pi_mix[k] + PI_FLOOR);
                if (sk2 == 0.0) {
                    log_g[k] = log_pi_k - 0.5 * z_sq;
                } else {
                    double s_k = u_bar * sk2 + 1.0;
                    log_g[k] = log_pi_k - 0.5 * (log(s_k) + z_sq / s_k);
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

        double rel_dsig = fabs(sg->sigma2 - old_sigma2) / (old_sigma2 + 1e-15);
        if (!silent) {
            fprintf(stderr, "var_fixed iter %3d: max_delta=%.2e  sigma2=%.4f"
                    "  h2_fixed=%.4f  rel_dsig=%.2e  pi=[",
                    iter + 1, max_delta, sg->sigma2, sg->h2_fixed, rel_dsig);
            for (int k = 0; k < K; k++)
                fprintf(stderr, " %.3f", sg->pi_mix[k]);
            fprintf(stderr, " ]\n");
        }

        old_sigma2 = sg->sigma2;
        if (rel_dsig < VAR_FIXED_TOL) break;
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
    const double        *ws2;        /* E[w_s^2] for each pruned segment, length S */
    double               sigma2;
    int                  n;
    int                  S;          /* pruned segment count (length of ws2) */
    int                  S_orig;     /* original tree count (for tskit iteration) */
    const int           *prune_map;  /* length S_orig: orig tree → pruned seg */
    const int           *hap_to_indiv;
    double              *work;      /* 2 * n_nodes scratch (borrowed) */
    tsk_id_t            *order_buf; /* n_nodes scratch (borrowed) */
    double              *gv;        /* n scratch (borrowed) */
} v_matvec_data_t;

static void v_matvec(double *out, const double *in, int n, void *data)
{
    v_matvec_data_t *d = data;

    for (int i = 0; i < n; i++) out[i] = d->sigma2 * in[i];

    /* Skip tree iteration entirely when all segment weights are zero */
    int has_random = 0;
    for (int s = 0; s < d->S; s++) {
        if (d->ws2[s] != 0.0) { has_random = 1; break; }
    }
    if (!has_random) return;

    tsk_tree_t tree;
    int ret = tsk_tree_init(&tree, d->ts, 0);
    if (ret < 0) {
        fprintf(stderr, "v_matvec: tsk_tree_init: %s\n", tsk_strerror(ret));
        return;
    }
    int s = 0;
    for (ret = tsk_tree_first(&tree); ret == TSK_TREE_OK && s < d->S_orig;
         ret = tsk_tree_next(&tree), s++) {
        int sp = d->prune_map[s];
        if (d->ws2[sp] == 0.0) continue;
        tree_cov_matvec(d->gv, in, d->n, &tree, d->work, d->order_buf, d->hap_to_indiv);
        for (int i = 0; i < n; i++) out[i] += d->ws2[sp] * d->gv[i];
    }
    if (ret < 0)
        fprintf(stderr, "v_matvec: tree error: %s\n", tsk_strerror(ret));
    tsk_tree_free(&tree);
}

/* Compute A[sp] = r_tilde^T G_sp r_tilde for all S pruned segments.
 * Iterates S_orig original trees; only the representative (first) tree of
 * each pruned segment is evaluated. */
static void compute_quad_all(const tsk_treeseq_t *ts, const int *hap_to_indiv,
                              int n, int S, int S_orig, const int *prune_map,
                              const double *r_tilde,
                              double *A,
                              double *work, tsk_id_t *order_buf, double *gv)
{
    for (int sp = 0; sp < S; sp++) A[sp] = 0.0;

    tsk_tree_t tree;
    int ret = tsk_tree_init(&tree, ts, 0);
    if (ret < 0) {
        fprintf(stderr, "compute_quad_all: tsk_tree_init: %s\n", tsk_strerror(ret));
        return;
    }
    int s = 0;
    for (ret = tsk_tree_first(&tree); ret == TSK_TREE_OK && s < S_orig;
         ret = tsk_tree_next(&tree), s++) {
        int sp = prune_map[s];
        tree_cov_matvec(gv, r_tilde, n, &tree, work, order_buf, hap_to_indiv);
        double a = 0.0;
        for (int i = 0; i < n; i++) a += r_tilde[i] * gv[i];
        A[sp] += a;   /* accumulate all trees in bin */
    }
    if (ret < 0)
        fprintf(stderr, "compute_quad_all: tree error: %s\n", tsk_strerror(ret));
    tsk_tree_free(&tree);
}


void var_random(sargent_t *sg, int silent)
{
    int n      = sg->n;
    int S      = sg->S;
    int nprobe = sg->nprobe;

    /* Allocate ws2 if needed; initialise to tau2_mix[M-1]/2 so that the
     * first CG solve sees a non-trivial V and the REML score has the right
     * sign (positive signal, not dominated by the initialisation bias). */
    if (!sg->ws2) {
        sg->ws2 = malloc((size_t)S * sizeof(double));
        if (!sg->ws2) { fprintf(stderr, "var_random: out of memory\n"); return; }
        double ws2_init = sg->tau2_mix[sg->M - 1] * 0.5;
        for (int s = 0; s < S; s++) sg->ws2[s] = ws2_init;
    }

    /* Working arrays */
    double   *r        = malloc((size_t)n * sizeof(double));
    double   *r_tilde  = calloc((size_t)n, sizeof(double));
    double   *A        = malloc((size_t)S * sizeof(double));
    double   *B        = malloc((size_t)S * sizeof(double));

    tsk_size_t n_nodes   = tsk_treeseq_get_num_nodes(sg->ts);
    double    *work      = malloc(2 * n_nodes * sizeof(double));
    tsk_id_t  *order_buf = malloc(n_nodes * sizeof(tsk_id_t));
    double    *gv        = malloc((size_t)n * sizeof(double));

    cg_state_t *cg = cg_state_new(n, NULL);

    if (!r || !r_tilde || !A || !B || !work || !order_buf || !gv || !cg) {
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

    v_matvec_data_t vd = {
        .ts = sg->ts, .ws2 = sg->ws2, .sigma2 = sg->sigma2,
        .n = n, .S = S, .S_orig = sg->S_orig, .prune_map = sg->prune_map,
        .hap_to_indiv = sg->hap_to_indiv,
        .work = work, .order_buf = order_buf, .gv = gv,
    };
    hutch_opts_t hopts = {
        .nprobe = nprobe, .batch_size = 0, .use_gaussian = 0, .seed = 0,
    };

    if (!silent)
        fprintf(stderr, "var_random: using CG-Hutchinson path (n=%d)\n", n);

    double old_sigma2 = sg->sigma2;

    for (int iter = 0; iter < VAR_RANDOM_MAX_ITER; iter++) {

        double A_0, B_0 = 0.0;

        vd.sigma2 = sg->sigma2;

        cg_solve(cg, v_matvec, &vd, r, r_tilde);

        A_0 = 0.0;
        for (int i = 0; i < n; i++) A_0 += r_tilde[i] * r_tilde[i];

        compute_quad_all(sg->ts, sg->hap_to_indiv, n, S, sg->S_orig, sg->prune_map,
                         r_tilde, A, work, order_buf, gv);

        hutch_trace_all(cg, v_matvec, &vd, sg->ts, sg->hap_to_indiv,
                        n, S, sg->S_orig, sg->prune_map, B, &B_0, &hopts);

        /* --- REML EM update for ws2 ----------------------------------------
         * Per-segment multiplicative update: ws2[s]_new = ws2[s]_old * A[s] / B[s].
         * Fixed point: A[s] = B[s] for all s, i.e. each segment's REML score = 0.
         * ws2 is floored at ws2_floor to prevent collapse; the Hutchinson noise
         * (~32% per segment with nprobe=10) will keep ws2[s] from reaching exact 0. */
        double ws2_max   = var_y / fmax(S, 1);
        double ws2_floor = 1e-12 * ws2_max;
        double sum_A = 0.0, sum_B = 0.0, ws2_mean_old = 0.0, ws2_mean_new = 0.0;
        for (int s = 0; s < S; s++) ws2_mean_old += sg->ws2[s];
        ws2_mean_old /= fmax(S, 1);

        for (int s = 0; s < S; s++) {
            sum_A += A[s];
            sum_B += B[s];
            if (B[s] > 1e-30)
                sg->ws2[s] *= A[s] / B[s];
            if (sg->ws2[s] < ws2_floor) sg->ws2[s] = ws2_floor;
            if (sg->ws2[s] > ws2_max)   sg->ws2[s] = ws2_max;
            ws2_mean_new += sg->ws2[s];
        }
        ws2_mean_new /= fmax(S, 1);

        double max_rel_ws2 = fabs(ws2_mean_new - ws2_mean_old) /
                             (ws2_mean_old + 1e-30);

        /* --- Update sigma2 (REML EM step) ---------------------------------- */
        if (B_0 > 0.0) {
            sg->sigma2 *= A_0 / B_0;
            if (sg->sigma2 < 1e-10 * var_y) sg->sigma2 = 1e-10 * var_y;
            if (sg->sigma2 > var_y)         sg->sigma2 = var_y;
        }

        sg->h2_random = (var_y > 0.0) ? fmax(0.0, 1.0 - sg->sigma2/var_y - sg->h2_fixed) : 0.0;
        sg->h2        = sg->h2_fixed + sg->h2_random;

        double rel_change = fabs(sg->sigma2 - old_sigma2) / (old_sigma2 + 1e-15);

        if (!silent) {
            fprintf(stderr, "var_random iter %3d: sigma2=%.4f  h2_rand=%.4f"
                    "  ws2_mean=%.3e  sumA=%.3e  sumB=%.3e"
                    "  B0=%.2e  rel_dsig=%.2e  rel_ws2=%.2e\n",
                    iter + 1, sg->sigma2, sg->h2_random,
                    ws2_mean_new, sum_A, sum_B,
                    B_0, rel_change, max_rel_ws2);
        }

        old_sigma2 = sg->sigma2;
        /* Convergence: sigma2 stable for at least 2 consecutive iterations.
         * ws2 oscillates at the Hutchinson noise floor and is not checked here. */
        if (iter >= 1 && rel_change < VAR_RANDOM_TOL)
            break;
    }

cleanup:
    cg_state_free(cg);
    free(r);
    free(r_tilde);
    free(A);
    free(B);
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

        double u_bar_post = (double)n / sg->sigma2;
        double z_sq = t_j * t_j / u_bar_post;
        double lmax = -1e300;
        for (int k = 0; k < K; k++) {
            double sk2 = sg->sigma2_mix[k];
            double log_pi_k = log(sg->pi_mix[k] + PI_FLOOR);
            if (sk2 == 0.0) {
                log_g[k] = log_pi_k - 0.5 * z_sq;
            } else {
                double s_k = u_bar_post * sk2 + 1.0;
                log_g[k] = log_pi_k - 0.5 * (log(s_k) + z_sq / s_k);
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
