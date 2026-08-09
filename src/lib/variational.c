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
#include "gall.h"

/* Default BayesR mixture grid: null + three effect-size classes.
 * sigma2_mix values are set relative to var(y) in sargent_new. */
#define VAR_DEFAULT_K        4
#define VAR_DEFAULT_M        2       /* spike-and-slab for random effects: null + slab */
#define VAR_FIXED_MAX_ITER   500
#define VAR_FIXED_POLYAK_BURN 20     /* sweeps skipped before Polyak-averaging beta (opt-in) */
#define VAR_RANDOM_MAX_ITER       100  /* total EM iterations (burn-in + averaging), cold-start */
#define VAR_RANDOM_MAX_ITER_WARM   50  /* total EM iterations when EMMA-warm-started */
#define FUSED_LASSO_PROX_MAX  200    /* coord-descent sweeps per prox call */
#define FUSED_LASSO_PROX_TOL  1e-6   /* sweep-max delta below which prox converges */
/* Trust-region cap on the Stage-1 fixed-effect variance: var(Xβ) ≤ this × var_y
 * each sweep. On mixed-ancestry training the BayesR coordinate descent otherwise
 * fits allele-frequency-divergent SNPs to separate the populations (confounding-
 * by-ancestry); the betas escape the prior (max|β| ~40× normal), σ² inflates to
 * ~4× var_y, and h²_snp = 1 − σ²/var_y diverges to ≈ −3. Capping var(Xβ) rescales
 * the betas back onto the feasible set; it only binds on the pathological
 * (overfit) cells — well-conditioned cells (var(Xβ) < cap·var_y) are untouched —
 * and recovers the degenerate cells to a sane h²_snp near the simulated truth.
 * This is a pragmatic regulariser, NOT a root-cause fix: the principled fix is
 * ancestry covariates in Stage 1 (see issue #22 / ROADMAP). Env override:
 * SARGENT_S2_VFCAP (0 disables). */
#define VAR_FIXED_VFCAP      0.6
/* Relaxed trust region used when ELBO-selected multi-start is active; see the
 * rationale at its use site in var_fixed_run. A 0.9 cap never binds on the
 * 500 kb h2=0.8 grid (identical to no cap) but still catches genuine blow-up. */
#define VAR_FIXED_VFCAP_MULTISTART 0.9

/* Block 2.6 — fused-lasso proximal helpers (sparsity-on-levels +
 * sparsity-on-differences on theta = log ws² − mean log ws²).
 *
 * The penalty is
 *     lambda1 * sum_s |theta_s|
 *   + lambda2 * sum_{s>0} |theta_s − theta_{s-1}|
 * (positional fusion weights; ARG-aware weights deferred to 2.6.x).
 *
 * The proximal step applied after each var_random EM update is:
 *     theta ← argmin { 0.5 ||theta − theta_in||² + lambda1·||theta||_1 +
 *                       lambda2·TV(theta) }
 * We solve it via sequential application of (TV-prox by cyclic
 * coordinate descent) then (soft-threshold for L1). This isn't the exact
 * joint ADMM solution, but the gap is small at typical λ scales and the
 * v1 implementation can be upgraded if 2.6 acceptance tests demand it. */

static inline double soft_threshold(double x, double lambda)
{
    if (x >  lambda) return x - lambda;
    if (x < -lambda) return x + lambda;
    return 0.0;
}

/* Single-coordinate TV update: minimise
 *     0.5 (x − y)² + lambda·(|x − a| + |x − b|)
 * with respect to x, given fixed neighbours a, b. */
static inline double tv_update_internal(double y, double a, double b, double lambda)
{
    double lo = a < b ? a : b;
    double hi = a < b ? b : a;
    if (y <  lo - 2.0 * lambda) return y + 2.0 * lambda;
    if (y <  lo)                return lo;
    if (y <= hi)                return y;
    if (y <= hi + 2.0 * lambda) return hi;
    return y - 2.0 * lambda;
}

/* Single-coordinate TV update at a boundary (only one neighbour a):
 *     min 0.5 (x − y)² + lambda·|x − a| */
static inline double tv_update_boundary(double y, double a, double lambda)
{
    if (y - a >  lambda) return y - lambda;
    if (y - a < -lambda) return y + lambda;
    return a;
}

/* 1D TV-L1 prox on theta (in-place). theta_in is the prox target;
 * theta_out receives the solution.  Cyclic coordinate descent with
 * convergence test, capped at FUSED_LASSO_PROX_MAX sweeps. */
static void fused_lasso_prox(const double *theta_in, double *theta_out,
                             int S, double lambda1, double lambda2)
{
    for (int s = 0; s < S; s++) theta_out[s] = theta_in[s];
    if (S < 1) return;

    if (lambda2 > 0.0 && S > 1) {
        for (int it = 0; it < FUSED_LASSO_PROX_MAX; it++) {
            double max_delta = 0.0;
            for (int s = 0; s < S; s++) {
                double y = theta_in[s];
                double prev = theta_out[s];
                double nx;
                if (s == 0) {
                    nx = tv_update_boundary(y, theta_out[1], lambda2);
                } else if (s == S - 1) {
                    nx = tv_update_boundary(y, theta_out[S - 2], lambda2);
                } else {
                    nx = tv_update_internal(y, theta_out[s - 1],
                                            theta_out[s + 1], lambda2);
                }
                theta_out[s] = nx;
                double d = fabs(nx - prev);
                if (d > max_delta) max_delta = d;
            }
            if (max_delta < FUSED_LASSO_PROX_TOL) break;
        }
    }

    if (lambda1 > 0.0) {
        for (int s = 0; s < S; s++)
            theta_out[s] = soft_threshold(theta_out[s], lambda1);
    }
}
/* Iteration cap for the Stage-2 V-solves (probe + main CG). The ARG covariance
 * V = sigma2 I + sum_s ws2[s] G_s is ill-conditioned for diverged-ancestry
 * cells (population-structure block of clustered small eigenvalues; tree
 * kinship has a CONSTANT diagonal, so Jacobi preconditioning gives no help).
 * Uncapped (limit = n), each probe CG solve runs the full n iterations on such
 * cells and one EM iteration takes ~100 s, so Stage 2 ran ~4 h (and got killed)
 * on cells whose well-conditioned siblings finished in ~2.5 h. The Hutchinson
 * trace estimator tolerates inexact solves (it already uses probe tol 1e-3, and
 * its statistical error is O(1/sqrt(nprobe))); capping at this value reproduces
 * the uncapped h2_arg / h2_snp / sigma2 / ws2 to 5 significant figures on
 * validation cells while bounding per-iteration cost. Overridable per run via
 * the SARGENT_CG_MAXIT environment variable. */
#define VAR_RANDOM_CG_MAXIT   64
#define VAR_RANDOM_BURN_IN        40  /* iterations before Polyak-Ruppert averaging, cold-start */
#define VAR_RANDOM_BURN_IN_WARM   15  /* burn-in when EMMA-warm-started */
#define VAR_FIXED_TOL        1e-4
#define PI_FLOOR             1e-10   /* lower bound on pi/phi_mix to avoid log(0) */
#define SCORE_CAP            50.0    /* clip tau^2*(A-B)/2 before softmax */
#define VAR_N_SEGMENTS       0       /* 0 = no pruning; use all original trees */

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
    sg->tau2_mix[0] = 0.0;           /* null component */
    sg->tau2_mix[1] = 1.0 * seg_scale;  /* slab component */

    /* Uniform initial responsibilities and proportions */
    for (int m = 0; m < sg->M; m++)
        sg->phi_mix[m] = 1.0 / sg->M;
    for (int s = 0; s < sg->S; s++)
        for (int m = 0; m < sg->M; m++)
            sg->rho[s * sg->M + m] = 1.0 / sg->M;

    sg->nprobe     = 10;   /* default; override via sg->nprobe before calling var_random */
    sg->phi_alpha  = 1.0;  /* uniform Beta(1,1) hyperprior on phi */
    sg->phi_beta   = 1.0;

    /* Alternative-mode flags; off by default */
    sg->uniform_arg = 0;
    sg->em_iters       = 0;
    sg->fused_lasso_lambda1 = 0.0;
    sg->fused_lasso_lambda2 = 0.0;

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
/* xorshift64* — small deterministic RNG for restart sweep-order permutations.
 * Self-contained so restarts reproduce from their seed alone and do not
 * disturb any global rand() state. */
static uint64_t vf_rng_next(uint64_t *s)
{
    uint64_t x = *s;
    x ^= x >> 12; x ^= x << 25; x ^= x >> 27;
    *s = x;
    return x * 2685821657736338717ULL;
}

static void vf_shuffle(int *a, int m, uint64_t *s)
{
    for (int i = m - 1; i > 0; i--) {
        int j = (int)(vf_rng_next(s) % (uint64_t)(i + 1));
        int t = a[i]; a[i] = a[j]; a[j] = t;
    }
}

/* Evaluate the mean-field variational objective (ELBO) at the current
 * variational parameters.  Used only to RANK restarts; the identical criterion
 * is applied to each, so the ranking is meaningful even though the
 * responsibilities are recomputed here from the committed beta rather than
 * carried out of the sweep (they differ slightly under Polyak averaging).
 *
 *   ELBO = -n/2 log(2 pi s2) - (||r||^2 + n * sum_j Var_q[beta_j]) / (2 s2)
 *          + sum_j sum_k    gamma_jk log(pi_k / gamma_jk)
 *          + sum_j sum_k>0  gamma_jk [ 1/2 log(v_jk/s_k^2) + 1/2
 *                                      - (mu_jk^2 + v_jk)/(2 s_k^2) ]
 *
 * with standardised genotypes (x_j'x_j = n).  Two genotype passes.
 * Returns -HUGE_VAL on allocation failure. */
static double var_fixed_elbo(sargent_t *sg, double *g, double *r)
{
    const int n = sg->n, p = sg->p, K = sg->K;
    const double se = sg->sigma2;
    if (!(se > 0.0)) return -HUGE_VAL;
    const double u_bar = (double)n / se;

    /* pass 1: residual from the committed beta */
    memcpy(r, sg->y, (size_t)n * sizeof(double));
    geno_rewind(sg->gr);
    for (int j = 0; j < p; j++) {
        if (geno_read_variant(sg->gr, g, GENO_STANDARDIZE) <= 0) break;
        double b = sg->beta[j];
        if (b != 0.0) for (int i = 0; i < n; i++) r[i] -= b * g[i];
    }

    double *lg = malloc((size_t)K * sizeof(double));
    double *gm = malloc((size_t)K * sizeof(double));
    if (!lg || !gm) { free(lg); free(gm); return -HUGE_VAL; }

    double var_sum = 0.0, kl = 0.0;

    /* pass 2: per-SNP variational terms at the final hyperparameters */
    geno_rewind(sg->gr);
    for (int j = 0; j < p; j++) {
        if (geno_read_variant(sg->gr, g, GENO_STANDARDIZE) <= 0) break;
        double bj = sg->beta[j];
        double xr = 0.0;
        for (int i = 0; i < n; i++) xr += g[i] * (r[i] + bj * g[i]);
        double t_j  = xr / se;
        double z_sq = t_j * t_j / u_bar;

        double lmax = -1e300;
        for (int k = 0; k < K; k++) {
            double sk2 = sg->sigma2_mix[k];
            double lpi = log(sg->pi_mix[k] + PI_FLOOR);
            if (sk2 == 0.0) lg[k] = lpi - 0.5 * z_sq;
            else {
                double s_k = u_bar * sk2 + 1.0;
                lg[k] = lpi - 0.5 * (log(s_k) + z_sq / s_k);
            }
            if (lg[k] > lmax) lmax = lg[k];
        }
        double sum_g = 0.0;
        for (int k = 0; k < K; k++) { gm[k] = exp(lg[k] - lmax); sum_g += gm[k]; }
        for (int k = 0; k < K; k++) gm[k] /= sum_g;

        double mean = 0.0, m2 = 0.0;
        for (int k = 0; k < K; k++) {
            if (gm[k] > 0.0)
                kl += gm[k] * (log(sg->pi_mix[k] + PI_FLOOR) - log(gm[k] + 1e-300));
            double sk2 = sg->sigma2_mix[k];
            if (sk2 == 0.0) continue;
            double v  = 1.0 / (u_bar + 1.0 / sk2);
            double mu = t_j * v;
            mean += gm[k] * mu;
            m2   += gm[k] * (mu * mu + v);
            kl   += gm[k] * (0.5 * log(v / sk2) + 0.5 - (mu * mu + v) / (2.0 * sk2));
        }
        var_sum += m2 - mean * mean;
    }
    free(lg); free(gm);

    double rr = 0.0;
    for (int i = 0; i < n; i++) rr += r[i] * r[i];
    double ll = -0.5 * (double)n * log(2.0 * M_PI * se)
                - (rr + (double)n * var_sum) / (2.0 * se);
    return ll + kl;
}

/* One Stage-1 fit.  randomized=0 keeps the historical sequential sweep order
 * (and the streaming genotype read) so a single non-randomized run is
 * byte-identical to <=0.4.0.  randomized=1 permutes the sweep order each
 * sweep, which is what lets restarts explore different basins. */
static void var_fixed_run(sargent_t *sg, int silent, uint64_t rseed, int randomized,
                          int multistart)
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
    int    *order     = randomized ? malloc((size_t)p * sizeof(int)) : NULL;
    if (!g || !r || !log_g || !gamma || !gamma_sum || (randomized && !order)) {
        fprintf(stderr, "var_fixed: out of memory\n");
        free(g); free(r); free(log_g); free(gamma); free(gamma_sum); free(order);
        return;
    }
    uint64_t rng_state = rseed ? rseed : 0x9E3779B97F4A7C15ULL;
    if (randomized) for (int j = 0; j < p; j++) order[j] = j;

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

    /* Fixed-effect variance cap (see VAR_FIXED_VFCAP above). Env override;
     * SARGENT_S2_VFCAP=0 disables. */
    /* The 0.6 clamp is tight enough to BIND on h2=0.8 traits, where it
     * under-scales beta and costs ~0.015 held-out r2 on every well-behaved
     * cell.  It earns that cost only because the single-start optimiser can
     * collapse into a dead basin (observed on 1 of 10 WGS h08 seeds:
     * h2_snp 0.80 -> 0.20, r2 0.74 -> 0.005 with the cap removed).
     * ELBO-selected multi-start removes that failure mode at the source, so
     * with restarts on the cap relaxes to a pure divergence guard that never
     * binds in practice.  Measured (10 seeds, dens3800, h2=0.8): restarts
     * with a 0.9 cap gain +0.0203 held-out r2 over the 0.6 default, 10/10
     * seeds, min 0.7173 (no collapse); a 0.9 cap and no cap are identical.
     * At h2=0.5 the cap never binds and the change is null (+0.0022). */
    double vf_cap = multistart ? VAR_FIXED_VFCAP_MULTISTART : VAR_FIXED_VFCAP;
    { const char *e = getenv("SARGENT_S2_VFCAP"); if (e) { double v = atof(e); if (v >= 0.0) vf_cap = v; } }

    /* Learn the BayesR mixture scale (EM on sigma2_mix each sweep, keeping the
     * relative grid fixed). OFF by default so existing results are byte-identical;
     * SARGENT_LEARN_SCALE=1 enables it. The frozen sigma2_mix set in sargent_new
     * ([0,1e-4,1e-3,1e-2]*var_y) does not adapt to the trait's effect scale and so
     * over-shrinks large effects (~8x below the data-driven scale on WGS cells,
     * per the sargent-paper D1 analysis). NOTE (D1 follow-up): learning the scale
     * turned out to be r2-NEUTRAL on the real estimator (held-out r2 gain
     * -0.0014, 95% CI [-0.0051,+0.0023]); it is a stability/h2_snp-calibration
     * option (repairs escaping-beta degeneracy on cross-ancestry cells), not a
     * prediction fix. Kept OFF by default for that reason. */
    int learn_scale = 0;
    { const char *e = getenv("SARGENT_LEARN_SCALE"); if (e && atoi(e) != 0) learn_scale = 1; }

    /* Polyak-Ruppert averaging of beta over post-burn-in sweeps. var_fixed's
     * stopping rule watches sigma2 (rel_dsig), but at p >> n beta can settle into a
     * limit cycle (mixture components flip-flopping) whose amplitude never decays
     * while sigma2 converges — so the run halts on a criterion blind to it (observed:
     * max_delta plateaus ~23% of max|beta| through exit). Averaging the tail iterate
     * is the standard fix (var_random already does this for Stage 2). ON by default
     * as of 0.4.0, validated on 10 WGS arg0 h08 seeds: the averaged beta leaves a
     * smaller training residual than the final iterate on 10/10 seeds (h2_snp
     * +0.0059, 95% CI [+0.0042,+0.0076], p=2e-5), confirming the exit point is an
     * off-centre sample of the cycle. The held-out r2 effect is small and marginal
     * (+0.0019, n=10, p=0.047), so this is an estimator-correctness fix rather than
     * a prediction fix. SARGENT_POLYAK_BETA=0 restores the pre-0.4.0 final-iterate
     * behaviour (byte-identical to <=0.3.0). */
    int polyak = 1;
    { const char *e = getenv("SARGENT_POLYAK_BETA"); if (e) polyak = (atoi(e) != 0); }
    double *beta_sum = NULL; long navg = 0;
    if (polyak) { beta_sum = calloc((size_t)p, sizeof(double)); if (!beta_sum) polyak = 0; }

    /* Convergence knobs (diagnostic; defaults reproduce <=0.4.0 exactly).
     * The stopping rule watches sigma2, which the file's own Polyak note
     * records as blind to a beta limit cycle -- so whether each restart is
     * actually AT an optimum is testable only by tightening these. */
    int    max_iter = VAR_FIXED_MAX_ITER;
    double s1_tol   = VAR_FIXED_TOL;
    { const char *e = getenv("SARGENT_S1_MAXITER"); if (e) { int v = atoi(e); if (v > 0) max_iter = v; } }
    { const char *e = getenv("SARGENT_S1_TOL");     if (e) { double v = atof(e); if (v > 0) s1_tol = v; } }

    for (int iter = 0; iter < max_iter; iter++) {

        /* u_bar = n/σ² is constant across all SNPs for standardised genotypes */
        double u_bar = (double)n / sg->sigma2;

        double max_delta = 0.0;
        memset(gamma_sum, 0, (size_t)K * sizeof(double));
        double scale_num = 0.0, scale_den = 0.0;   /* EM sufficient stats for the mixture scale */
        int n_read = 0;

        if (randomized) vf_shuffle(order, p, &rng_state);
        else            geno_rewind(sg->gr);

        for (int oi = 0; oi < p; oi++) {
            int j = randomized ? order[oi] : oi;
            if (randomized) geno_seek(sg->gr, j);
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
                /* EM stats for the mixture scale (non-null components only):
                 * factor = Σ γ_jk E[β_jk²]/σ_k² ÷ Σ γ_jk, applied after the sweep */
                scale_num    += gamma[k] * (mu_jk * mu_jk + v_jk) / sg->sigma2_mix[k];
                scale_den    += gamma[k];
            }

            /* Damped update: β ← ½β_old + ½β_posterior.
             * Prevents the BayesR mixture-snap cascade (same root cause as
             * var_joint Phase 4.5 fix, commit bdce0e6): a SNP flipping
             * components makes a discontinuous residual jump that destabilises
             * σ².  Half-stepping makes each transition gradual; convergence
             * takes more iterations but stays stable when p >> n.           */
            double beta_damped = 0.5 * beta_old + 0.5 * beta_new;
            sg->beta[j] = beta_damped;

            /* Update residual using the damped beta actually committed */
            for (int i = 0; i < n; i++) r[i] -= beta_damped * g[i];

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

        /* Learn the mixture scale (EM), preserving the relative grid: multiply
         * every non-null σ_k² by the same factor = Σγ E[β²]/σ_k² ÷ Σγ. Guards
         * keep it in the regime the freeze was protecting: a per-sweep factor
         * clamp [0.25,4] (no single-step snap) and an absolute cap σ_k² ≤ var_y
         * (no single-SNP effect exceeds the phenotypic variance); the VFCAP
         * trust region below is the final backstop on var(Xβ). */
        if (learn_scale && scale_den > 0.0) {
            double factor = scale_num / scale_den;
            if (factor < 0.25) factor = 0.25;
            if (factor > 4.0)  factor = 4.0;
            for (int k = 0; k < K; k++) {
                if (sg->sigma2_mix[k] <= 0.0) continue;
                double s2 = sg->sigma2_mix[k] * factor;
                sg->sigma2_mix[k] = (s2 > var_y) ? var_y : s2;
            }
        }

        /* Trust region on the fixed effect: the SNP-explained component
         * var(Xβ) cannot exceed the total variance for a valid fit; an excess
         * is the cross-ancestry overfit (betas escaping the prior). Rescale β
         * (and the residual) onto var(Xβ) ≤ vf_cap·var_y. Since Xβ = y − r,
         * the rescale needs no genotype pass. */
        if (vf_cap > 0.0) {
            double fit_ss = 0.0;
            for (int i = 0; i < n; i++) { double f = sg->y[i] - r[i]; fit_ss += f * f; }
            double vf = fit_ss / n, cap = vf_cap * var_y;
            if (vf > cap && vf > 0.0) {
                double s = sqrt(cap / vf);
                for (int j = 0; j < p; j++) sg->beta[j] *= s;
                for (int i = 0; i < n; i++) r[i] = (1.0 - s) * sg->y[i] + s * r[i];
            }
        }

        /* Update residual variance: σ² = r^T r / n, with damped upward moves.
         * Downward moves are accepted in full (they represent genuine descent).
         * Upward moves are half-stepped: σ²_new ← ½σ²_old + ½(rr/n).
         * Without this, a mixture-snap that zeroes several betas can double
         * σ² in one iteration, driving the cascade that Phase 4.5 (bdce0e6)
         * diagnosed for var_joint — the same instability affects var_fixed at
         * moderate n (observed: n=1000, p=3340 diverges without the cap).
         *
         * Note: var_fixed uses the MLE direct update ||r||²/n (bounded by
         * var(y)), so asymmetric damping suffices.  var_joint uses a Newton
         * step σ² *= A_0/B_0 whose ratio can be unbounded; it needs
         * additional symmetric + clip damping.  See ROADMAP "σ² damping
         * convention" subsection for the design rationale. */
        double rr = 0.0;
        for (int i = 0; i < n; i++) rr += r[i] * r[i];
        double sigma2_raw = rr / n;
        if (sigma2_raw < 1e-10 * var_y) sigma2_raw = 1e-10 * var_y;
        sg->sigma2 = (sigma2_raw < sg->sigma2)
                     ? sigma2_raw                             /* downward: full step */
                     : 0.5 * sg->sigma2 + 0.5 * sigma2_raw; /* upward:   half step */

        /* Update h2_snp: fraction of phenotypic variance explained by X β̂ */
        sg->h2_snp  = (var_y > 0.0) ? 1.0 - sg->sigma2 / var_y : 0.0;

        double rel_dsig = fabs(sg->sigma2 - old_sigma2) / (old_sigma2 + 1e-15);
        if (!silent) {
            fprintf(stderr, "var_fixed iter %3d: max_delta=%.2e  sigma2=%.4f"
                    "  h2_snp=%.4f  rel_dsig=%.2e  pi=[",
                    iter + 1, max_delta, sg->sigma2, sg->h2_snp, rel_dsig);
            for (int k = 0; k < K; k++)
                fprintf(stderr, " %.3f", sg->pi_mix[k]);
            fprintf(stderr, " ]\n");
        }

        if (polyak && iter >= VAR_FIXED_POLYAK_BURN) {
            for (int j = 0; j < p; j++) beta_sum[j] += sg->beta[j];
            navg++;
        }

        old_sigma2 = sg->sigma2;
        if (rel_dsig < s1_tol) break;
    }

    /* Commit the Polyak average and recompute sigma2/h2_snp from it for consistency. */
    if (polyak && navg > 0) {
        for (int j = 0; j < p; j++) sg->beta[j] = beta_sum[j] / (double)navg;
        memcpy(r, sg->y, (size_t)n * sizeof(double));
        geno_rewind(sg->gr);
        for (int j = 0; j < p; j++) {
            if (geno_read_variant(sg->gr, g, GENO_STANDARDIZE) <= 0) break;
            double b = sg->beta[j];
            if (b != 0.0) for (int i = 0; i < n; i++) r[i] -= b * g[i];
        }
        double rr = 0.0;
        for (int i = 0; i < n; i++) rr += r[i] * r[i];
        double s2 = rr / n;
        if (s2 < 1e-10 * var_y) s2 = 1e-10 * var_y;
        sg->sigma2 = s2;
        sg->h2_snp = (var_y > 0.0) ? 1.0 - sg->sigma2 / var_y : 0.0;
    }
    free(beta_sum);

    free(g);
    free(r);
    free(log_g);
    free(gamma);
    free(gamma_sum);
    free(order);
}

/* Stage 1 entry point: optional ELBO-selected multi-start.
 *
 * SARGENT_STAGE1_RESTARTS=N (default 1) runs N fits and keeps the one with the
 * highest ELBO.  N=1 takes the historical sequential-sweep path and is
 * byte-identical to <=0.4.0; N>1 randomises the sweep order per restart.
 *
 * Motivation (sargent-paper D1, 2026-08-09): the mean-field objective is
 * multimodal on high-coverage 500 kb cells.  A python mean-field mimic gained
 * +0.0353 held-out r2 (mean over 3 seeds) from best-of-6 ELBO-selected
 * restarts, and a single sequential fit collapsed outright on 1 of 10 seeds
 * (h2_snp 0.80 -> 0.20) once the var(X beta) trust region was relaxed.  Both
 * are basin-selection failures rather than a limitation of mean-field itself.
 */
/* ---------------------------------------------------------------------------
 * Stage 1, alternative inference: BayesR Gibbs (SSVS)
 *
 * Same model var_fixed fits variationally --
 *   y = X beta + e,  e ~ N(0, s2),  beta_j ~ sum_k pi_k N(0, s2_mix[k]),
 *   s2_mix[0] = 0 --
 * but sampled rather than approximated by a factorised posterior, so it does
 * not inherit mean-field's LD failure mode.  beta is reported as the
 * post-burn-in posterior mean.
 *
 * Why: on 500 kb WGS cells the variational fit sits ~0.017 held-out r2 below
 * BayesR (GCTB) even with ELBO-selected multi-start, and the sargent-paper D1
 * prototype put a Gibbs sampler at GCTB parity (+0.0010 over 3 seeds).  Every
 * mean-field-level lever was exhausted first: restart count saturates, an
 * oracle selector is worth only +0.0052, scale-learning is null, and tighter
 * convergence is actively harmful.
 *
 * Per sweep this costs about one var_fixed sweep, so the run cost is set by
 * the iteration count rather than by convergence.
 * ------------------------------------------------------------------------- */

static double vf_unif(uint64_t *s)
{
    return (double)(vf_rng_next(s) >> 11) * (1.0 / 9007199254740992.0);  /* [0,1) */
}

static double vf_norm(uint64_t *s)
{
    double u1, u2;
    do { u1 = vf_unif(s); } while (u1 <= 1e-300);
    u2 = vf_unif(s);
    return sqrt(-2.0 * log(u1)) * cos(2.0 * M_PI * u2);
}

/* Marsaglia-Tsang gamma(shape, 1) */
static double vf_gamma_rand(uint64_t *s, double shape)
{
    if (shape <= 0.0) return 0.0;
    if (shape < 1.0) {
        double u = vf_unif(s);
        if (u <= 0.0) u = 1e-300;
        return vf_gamma_rand(s, shape + 1.0) * pow(u, 1.0 / shape);
    }
    double d = shape - 1.0 / 3.0, c = 1.0 / sqrt(9.0 * d);
    for (;;) {
        double x = vf_norm(s);
        double v = 1.0 + c * x;
        if (v <= 0.0) continue;
        v = v * v * v;
        double u = vf_unif(s);
        if (u <= 0.0) u = 1e-300;
        if (u < 1.0 - 0.0331 * x * x * x * x) return d * v;
        if (log(u) < 0.5 * x * x + d * (1.0 - v + log(v))) return d * v;
    }
}

static void var_fixed_gibbs(sargent_t *sg, int silent, int niter, int burn,
                            uint64_t seed, int sample_scale)
{
    const int n = sg->n, p = sg->p, K = sg->K;

    if (!sg->beta) {
        sg->beta = calloc((size_t)p, sizeof(double));
        if (!sg->beta) { fprintf(stderr, "var_fixed_gibbs: out of memory\n"); return; }
    }

    double *g        = malloc((size_t)n * sizeof(double));
    double *r        = malloc((size_t)n * sizeof(double));
    double *lw       = malloc((size_t)K * sizeof(double));
    double *cnt      = malloc((size_t)K * sizeof(double));
    double *relgrid  = malloc((size_t)K * sizeof(double));
    double *beta_sum = calloc((size_t)p, sizeof(double));
    double *pi_sum   = calloc((size_t)K, sizeof(double));
    if (!g || !r || !lw || !cnt || !relgrid || !beta_sum || !pi_sum) {
        fprintf(stderr, "var_fixed_gibbs: out of memory\n");
        free(g); free(r); free(lw); free(cnt); free(relgrid); free(beta_sum); free(pi_sum);
        return;
    }

    /* Relative mixture grid, held fixed; only its common scale is sampled. */
    for (int k = 0; k < K; k++) relgrid[k] = sg->sigma2_mix[k];

    const double var_y = vec_var(sg->y, n);
    uint64_t st = seed ? seed : 0xD1B54A32D192ED03ULL;

    /* Common scale factor on the relative grid: sigma2_mix[k] = relgrid[k]*f.
     * f_max keeps the largest component at or below var_y (no single SNP
     * explains more than the phenotype); f_min stops total collapse. */
    double relmax = 0.0;
    for (int k = 0; k < K; k++) if (relgrid[k] > relmax) relmax = relgrid[k];
    const double f_max = (relmax > 0.0 && var_y > 0.0) ? var_y / relmax : 1e6;
    const double f_min = 1e-6;
    double scale_f = 1.0;

    /* r = y - X beta (beta is 0 on a cold start) */
    memcpy(r, sg->y, (size_t)n * sizeof(double));
    geno_rewind(sg->gr);
    for (int j = 0; j < p; j++) {
        if (geno_read_variant(sg->gr, g, GENO_STANDARDIZE) <= 0) break;
        double b = sg->beta[j];
        if (b != 0.0) for (int i = 0; i < n; i++) r[i] -= b * g[i];
    }

    long nsamp = 0;
    double sigma2_sum = 0.0;

    int resync = 0;
    { const char *e = getenv("SARGENT_GIBBS_RESYNC"); if (e && atoi(e) != 0) resync = 1; }

    for (int it = 0; it < niter; it++) {
        /* Diagnostic: rebuild r = y - X beta from scratch instead of trusting
         * the incremental +=/-= maintained across ~p*niter updates. */
        if (resync) {
            memcpy(r, sg->y, (size_t)n * sizeof(double));
            geno_rewind(sg->gr);
            for (int j = 0; j < p; j++) {
                if (geno_read_variant(sg->gr, g, GENO_STANDARDIZE) <= 0) break;
                double b = sg->beta[j];
                if (b != 0.0) for (int i = 0; i < n; i++) r[i] -= b * g[i];
            }
        }
        memset(cnt, 0, (size_t)K * sizeof(double));
        double ss_scaled = 0.0;      /* sum_j beta_j^2 / relgrid[k_j], non-null only */
        long   m_nonnull = 0;

        geno_rewind(sg->gr);
        for (int j = 0; j < p; j++) {
            if (geno_read_variant(sg->gr, g, GENO_STANDARDIZE) <= 0) break;

            double bold = sg->beta[j];
            if (bold != 0.0) for (int i = 0; i < n; i++) r[i] += bold * g[i];  /* r_{-j} */

            double xr = 0.0;
            for (int i = 0; i < n; i++) xr += g[i] * r[i];

            /* component log-weights; the exp(-r'r/2s2) factor is common and cancels */
            double lmax = -1e300;
            for (int k = 0; k < K; k++) {
                double sk2 = sg->sigma2_mix[k];
                double lpi = log(sg->pi_mix[k] + PI_FLOOR);
                if (sk2 <= 0.0) {
                    lw[k] = lpi;
                } else {
                    double v  = 1.0 / ((double)n / sg->sigma2 + 1.0 / sk2);
                    double mu = v * xr / sg->sigma2;
                    lw[k] = lpi + 0.5 * log(v / sk2) + mu * mu / (2.0 * v);
                }
                if (lw[k] > lmax) lmax = lw[k];
            }
            double tot = 0.0;
            for (int k = 0; k < K; k++) { lw[k] = exp(lw[k] - lmax); tot += lw[k]; }

            double u = vf_unif(&st) * tot, acc = 0.0;
            int ksel = K - 1;
            for (int k = 0; k < K; k++) { acc += lw[k]; if (u <= acc) { ksel = k; break; } }

            double bnew = 0.0;
            if (sg->sigma2_mix[ksel] > 0.0) {
                double v  = 1.0 / ((double)n / sg->sigma2 + 1.0 / sg->sigma2_mix[ksel]);
                double mu = v * xr / sg->sigma2;
                bnew = mu + sqrt(v) * vf_norm(&st);
                if (relgrid[ksel] > 0.0) { ss_scaled += bnew * bnew / relgrid[ksel]; m_nonnull++; }
            }
            sg->beta[j] = bnew;
            cnt[ksel] += 1.0;
            if (bnew != 0.0) for (int i = 0; i < n; i++) r[i] -= bnew * g[i];
        }

        /* pi ~ Dirichlet(1 + counts) */
        double psum = 0.0;
        for (int k = 0; k < K; k++) { sg->pi_mix[k] = vf_gamma_rand(&st, 1.0 + cnt[k]); psum += sg->pi_mix[k]; }
        if (psum > 0.0) for (int k = 0; k < K; k++) sg->pi_mix[k] /= psum;

        /* Common mixture scale f, with sigma2_mix[k] = relgrid[k] * f.
         *
         *   f ~ Inv-Gamma(a0 + m/2, b0 + ss_scaled/2)
         *
         * The prior must be PROPER.  Under a flat prior (a0=b0=0) a single
         * large beta landing in the smallest component contributes
         * beta^2/1e-4 to ss_scaled, f explodes, every component saturates at
         * var_y, the mixture degenerates to one diffuse component and sigma2
         * runs away (observed: sigma2 -> 133 against var_y ~ 1).  a0=2, b0=1
         * is a weak proper prior with mean f=1, i.e. centred on SARGENT's
         * frozen grid with ~4 degrees of freedom -- the same role GCTB's
         * scaled-inverse-chi-square prior on sigma2_beta plays. */
        if (sample_scale) {
            const double a0 = 2.0, b0 = 1.0;
            double gs = vf_gamma_rand(&st, a0 + 0.5 * (double)m_nonnull);
            if (gs > 0.0) {
                double f = (b0 + 0.5 * ss_scaled) / gs;
                /* Bound f itself, never the components individually.  Clamping
                 * each sigma2_mix[k] separately desynchronises sampling from
                 * the update: beta would be drawn with variance var_y but
                 * scored as beta^2/relgrid[k] (1e-4), inflating ss_scaled
                 * ~1e4-fold and pinning the chain at the cap forever. */
                if (f < f_min) f = f_min;
                if (f > f_max) f = f_max;
                scale_f = f;
                for (int k = 0; k < K; k++) sg->sigma2_mix[k] = relgrid[k] * scale_f;
            }
        }

        /* sigma2 ~ Inv-Gamma(n/2, r'r/2), flat prior */
        double rr = 0.0;
        for (int i = 0; i < n; i++) rr += r[i] * r[i];
        double gs2 = vf_gamma_rand(&st, 0.5 * (double)n);
        if (gs2 > 0.0) {
            double s2 = (0.5 * rr) / gs2;
            if (s2 < 1e-10 * var_y) s2 = 1e-10 * var_y;
            sg->sigma2 = s2;
        }

        if (it >= burn) {
            for (int j = 0; j < p; j++) beta_sum[j] += sg->beta[j];
            for (int k = 0; k < K; k++) pi_sum[k]   += sg->pi_mix[k];
            sigma2_sum += sg->sigma2;
            nsamp++;
        }
        if (!silent && (it % 50 == 0 || it == niter - 1)) {
            double bmax = 0.0, bss = 0.0;
            for (int j = 0; j < p; j++) { double a = fabs(sg->beta[j]); if (a > bmax) bmax = a; bss += a*a; }
            fprintf(stderr, "gibbs %4d/%d: s2=%.4g pi0=%.4f m=%ld f=%.4g max|b|=%.4g sum b^2=%.4g s2mix3=%.4g\n",
                    it + 1, niter, sg->sigma2, sg->pi_mix[0], m_nonnull, scale_f, bmax, bss,
                    sg->sigma2_mix[K-1]);
        }
    }

    /* Posterior means, then sigma2/h2_snp recomputed from the mean beta so the
     * reported summary is coherent with the beta actually written out. */
    if (nsamp > 0) {
        for (int j = 0; j < p; j++) sg->beta[j] = beta_sum[j] / (double)nsamp;
        for (int k = 0; k < K; k++) sg->pi_mix[k] = pi_sum[k] / (double)nsamp;

        memcpy(r, sg->y, (size_t)n * sizeof(double));
        geno_rewind(sg->gr);
        for (int j = 0; j < p; j++) {
            if (geno_read_variant(sg->gr, g, GENO_STANDARDIZE) <= 0) break;
            double b = sg->beta[j];
            if (b != 0.0) for (int i = 0; i < n; i++) r[i] -= b * g[i];
        }
        double rr = 0.0;
        for (int i = 0; i < n; i++) rr += r[i] * r[i];
        double s2 = rr / n;
        if (s2 < 1e-10 * var_y) s2 = 1e-10 * var_y;
        sg->sigma2 = s2;
        sg->h2_snp = (var_y > 0.0) ? 1.0 - sg->sigma2 / var_y : 0.0;
    }

    if (!silent)
        fprintf(stderr, "gibbs: %ld post-burn-in samples, h2_snp=%.4f sigma2=%.4f\n",
                nsamp, sg->h2_snp, sg->sigma2);

    free(g); free(r); free(lw); free(cnt); free(relgrid); free(beta_sum); free(pi_sum);
}

/* Set by --stage1-gibbs / --gibbs-iter / --gibbs-burn (see sargent.c). */
static int s1_gibbs_on = 0, s1_gibbs_iter = 3000, s1_gibbs_burn = 1000;

void var_fixed_set_gibbs(int on, int niter, int burn)
{
    s1_gibbs_on = on;
    if (niter > 0) s1_gibbs_iter = niter;
    if (burn  >= 0) s1_gibbs_burn = burn;
}

void var_fixed(sargent_t *sg, int silent)
{
    if (s1_gibbs_on) {
        int burn = s1_gibbs_burn < s1_gibbs_iter ? s1_gibbs_burn : s1_gibbs_iter / 2;
        int freeze = 0;
        { const char *e = getenv("SARGENT_GIBBS_FREEZE_SCALE"); if (e && atoi(e) != 0) freeze = 1; }
        uint64_t sd = 0xD1B54A32D192ED03ULL;
        { const char *e = getenv("SARGENT_GIBBS_SEED"); if (e) sd = (uint64_t)strtoull(e, NULL, 10); }

        /* Warm-start the chain from the variational fit.  sargent_new leaves
         * pi uniform (1/K), i.e. (K-1)/K of SNPs non-null; at p > n a sampler
         * started there saturates immediately and sigma2 runs away (observed:
         * 1111/1596 non-null on sweep 1, sigma2 -> 139 against var_y ~ 1).
         * The variational fit supplies a sparse pi and a sensible sigma2 in
         * ~3 s, which is also a standard way to cut MCMC burn-in.
         * SARGENT_GIBBS_COLDSTART=1 skips it (diagnostic). */
        int cold = 0;
        { const char *e = getenv("SARGENT_GIBBS_COLDSTART"); if (e && atoi(e) != 0) cold = 1; }
        if (!cold) {
            if (!silent) fprintf(stderr, "gibbs: warm-starting from the variational fit\n");
            var_fixed_run(sg, silent, 0, 0, 0);
        }
        var_fixed_gibbs(sg, silent, s1_gibbs_iter, burn, sd, !freeze);
        return;
    }

    int nrestart = 1;
    { const char *e = getenv("SARGENT_STAGE1_RESTARTS");
      if (e) { int v = atoi(e); if (v >= 1) nrestart = v; } }

    if (nrestart == 1) {          /* historical path, bit-for-bit */
        var_fixed_run(sg, silent, 0, 0, 0);
        return;
    }

    const int n = sg->n, p = sg->p, K = sg->K;

    /* Snapshot the entering state so every restart begins identically. */
    double *beta0 = sg->beta ? malloc((size_t)p * sizeof(double)) : NULL;
    double *pi0   = malloc((size_t)K * sizeof(double));
    double *s2m0  = malloc((size_t)K * sizeof(double));
    double *best_beta = malloc((size_t)p * sizeof(double));
    double *best_pi   = malloc((size_t)K * sizeof(double));
    double *best_s2m  = malloc((size_t)K * sizeof(double));
    double *gs        = malloc((size_t)n * sizeof(double));
    double *rs        = malloc((size_t)n * sizeof(double));
    if (!pi0 || !s2m0 || !best_beta || !best_pi || !best_s2m || !gs || !rs) {
        fprintf(stderr, "var_fixed: out of memory (restarts); falling back to 1 fit\n");
        free(beta0); free(pi0); free(s2m0);
        free(best_beta); free(best_pi); free(best_s2m); free(gs); free(rs);
        var_fixed_run(sg, silent, 0, 0, 0);
        return;
    }
    if (beta0) memcpy(beta0, sg->beta, (size_t)p * sizeof(double));
    memcpy(pi0,  sg->pi_mix,     (size_t)K * sizeof(double));
    memcpy(s2m0, sg->sigma2_mix, (size_t)K * sizeof(double));
    const double sigma2_0 = sg->sigma2;

    double best_elbo = -HUGE_VAL, best_sigma2 = sg->sigma2, best_h2 = sg->h2_snp;
    int best_k = -1;

    for (int k = 0; k < nrestart; k++) {
        /* restore the entering state */
        if (beta0 && sg->beta) memcpy(sg->beta, beta0, (size_t)p * sizeof(double));
        else if (sg->beta)     memset(sg->beta, 0, (size_t)p * sizeof(double));
        memcpy(sg->pi_mix,     pi0,  (size_t)K * sizeof(double));
        memcpy(sg->sigma2_mix, s2m0, (size_t)K * sizeof(double));
        sg->sigma2 = sigma2_0;

        /* restart 0 is the historical sequential sweep, so multi-start can
         * never do worse than the default fit by more than ELBO mis-ranking */
        var_fixed_run(sg, silent, 0x9E3779B97F4A7C15ULL + 0x1000ULL * (uint64_t)(k + 1),
                      k == 0 ? 0 : 1, 1);

        double e = var_fixed_elbo(sg, gs, rs);
        /* Diagnostic: SARGENT_STAGE1_DUMP_RESTARTS=<prefix> writes each
         * restart's beta so the ELBO ranking can be scored against held-out
         * r2 offline (is selection, or search, the bottleneck?). Debug only. */
        { const char *dp = getenv("SARGENT_STAGE1_DUMP_RESTARTS");
          if (dp && *dp) {
              char path[4096];
              snprintf(path, sizeof path, "%s_restart%d.beta.tsv", dp, k);
              FILE *fp = fopen(path, "w");
              if (fp) {
                  fprintf(fp, "# elbo %.10g h2_snp %.10g sigma2 %.10g\n",
                          e, sg->h2_snp, sg->sigma2);
                  for (int j = 0; j < p; j++) fprintf(fp, "%.10g\n", sg->beta[j]);
                  fclose(fp);
              }
          } }
        if (!silent)
            fprintf(stderr, "var_fixed restart %d/%d: elbo=%.4f h2_snp=%.4f sigma2=%.4f\n",
                    k + 1, nrestart, e, sg->h2_snp, sg->sigma2);
        if (e > best_elbo) {
            best_elbo = e; best_k = k;
            memcpy(best_beta, sg->beta,       (size_t)p * sizeof(double));
            memcpy(best_pi,   sg->pi_mix,     (size_t)K * sizeof(double));
            memcpy(best_s2m,  sg->sigma2_mix, (size_t)K * sizeof(double));
            best_sigma2 = sg->sigma2; best_h2 = sg->h2_snp;
        }
    }

    if (best_k >= 0) {
        memcpy(sg->beta,       best_beta, (size_t)p * sizeof(double));
        memcpy(sg->pi_mix,     best_pi,   (size_t)K * sizeof(double));
        memcpy(sg->sigma2_mix, best_s2m,  (size_t)K * sizeof(double));
        sg->sigma2 = best_sigma2; sg->h2_snp = best_h2;
        if (!silent)
            fprintf(stderr, "var_fixed: selected restart %d/%d (elbo=%.4f)\n",
                    best_k + 1, nrestart, best_elbo);
    }

    free(beta0); free(pi0); free(s2m0);
    free(best_beta); free(best_pi); free(best_s2m); free(gs); free(rs);
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
    double              *gv;         /* n scratch (borrowed) */
    tree_cov_matvec_inc_state_t *inc; /* borrowed; long-lived incremental state */
    int                  use_span;   /* span-weight the kernel? (uniform arms yes, learned no) */
} v_matvec_data_t;

/* --- v_matvec via tree_walker_run --- */

typedef struct {
    double       *out;
    const double *ws2;
    int           n;
} v_matvec_visit_ctx_t;

static int v_matvec_pre(int s_orig, int sp, void *ctx)
{
    (void)s_orig;
    const v_matvec_visit_ctx_t *vc = ctx;
    return vc->ws2[sp] != 0.0;
}

static void v_matvec_visit(const tree_walk_state_t *st, void *ctx)
{
    v_matvec_visit_ctx_t *vc = ctx;
    double w = vc->ws2[st->sp] * st->span;   /* span-weighted kernel: Sum_s ws2_s span_s G_s */
    for (int i = 0; i < vc->n; i++) vc->out[i] += w * st->gv[i];
}

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

    v_matvec_visit_ctx_t vc = { .out = out, .ws2 = d->ws2, .n = n };
    tree_walker_run_inc(d->inc, d->S_orig, d->prune_map, in, d->gv,
                        v_matvec_pre, v_matvec_visit, &vc, d->use_span);
}

/* --- compute_quad_all via tree_walker_run --- */

typedef struct {
    double       *A;
    const double *r_tilde;
    int           n;
} quad_visit_ctx_t;

static void quad_visit(const tree_walk_state_t *st, void *ctx)
{
    quad_visit_ctx_t *qc = ctx;
    double a = 0.0;
    for (int i = 0; i < qc->n; i++) a += qc->r_tilde[i] * st->gv[i];
    qc->A[st->sp] += st->span * a;   /* span-weighted: r_tilde^T (span_s G_s) r_tilde */
}

/* Compute A[sp] = r_tilde^T G_sp r_tilde for all S pruned segments. */
static void compute_quad_all(tree_cov_matvec_inc_state_t *inc,
                              int n, int S, int S_orig, const int *prune_map,
                              const double *r_tilde,
                              double *A, double *gv, int use_span)
{
    for (int sp = 0; sp < S; sp++) A[sp] = 0.0;
    quad_visit_ctx_t qc = { .A = A, .r_tilde = r_tilde, .n = n };
    tree_walker_run_inc(inc, S_orig, prune_map, r_tilde, gv,
                        NULL, quad_visit, &qc, use_span);
}


void var_random(sargent_t *sg, int silent)
{
    int n      = sg->n;
    int S      = sg->S;
    int nprobe = sg->nprobe;

    /* Detect EMMA warm-start: ws2 pre-filled by the caller (e.g. via a prior
     * EMMA call in sargent.c) → tighter iteration bounds since we start near
     * convergence.  Cold-start: allocate and set a default initial value. */
    int was_warmstarted = (sg->ws2 != NULL);
    if (!sg->ws2) {
        sg->ws2 = malloc((size_t)S * sizeof(double));
        if (!sg->ws2) { fprintf(stderr, "var_random: out of memory\n"); return; }
        double ws2_init = sg->tau2_mix[sg->M - 1] * 0.5;
        for (int s = 0; s < S; s++) sg->ws2[s] = ws2_init;
    }
    int max_iter = was_warmstarted ? VAR_RANDOM_MAX_ITER_WARM : VAR_RANDOM_MAX_ITER;
    int burn_in  = was_warmstarted ? VAR_RANDOM_BURN_IN_WARM  : VAR_RANDOM_BURN_IN;

    /* Working arrays */
    double   *r          = malloc((size_t)n * sizeof(double));
    double   *r_tilde    = calloc((size_t)n,  sizeof(double));
    double   *A          = malloc((size_t)S * sizeof(double));
    double   *B          = malloc((size_t)S * sizeof(double));
    double   *ws2_accum  = calloc((size_t)S,  sizeof(double));

    /* Block 2.6 — fused-lasso prox scratch (allocated only when needed). */
    int use_fl = (sg->fused_lasso_lambda1 > 0.0 || sg->fused_lasso_lambda2 > 0.0);
    double   *fl_log_ws2   = use_fl ? malloc((size_t)S * sizeof(double)) : NULL;
    double   *fl_theta_in  = use_fl ? malloc((size_t)S * sizeof(double)) : NULL;
    double   *fl_theta_out = use_fl ? malloc((size_t)S * sizeof(double)) : NULL;

    double    *gv        = malloc((size_t)n * sizeof(double));

    /* Bound the Stage-2 CG work (see VAR_RANDOM_CG_MAXIT above). Cap at the
     * smaller of n and the default; SARGENT_CG_MAXIT overrides for tuning. */
    int cg_cap = (n < VAR_RANDOM_CG_MAXIT) ? n : VAR_RANDOM_CG_MAXIT;
    { const char *e = getenv("SARGENT_CG_MAXIT"); if (e && atoi(e) > 0) cg_cap = atoi(e); }
    cg_opts_t main_cg_opts = {
        .max_iter     = cg_cap,
        .tol          = 1e-6,
        .abs_tol      = 0.0,
        .diag_precond = NULL,
        .precond_fn   = NULL,
        .precond_data = NULL,
    };
    cg_state_t *cg = cg_state_new(n, &main_cg_opts);

    /* Wire G_all preconditioner for per-tree mode when precomputed.
     * Reduces κ(P⁻¹V) from O(√κ)≈64 to ~O(1), cutting CG iterations ~64→5.
     * Pinv is rebuilt each iteration (O(n)) from the current ws2_mean/sigma2. */
    int use_gall_precond = (sg->gall != NULL && !sg->uniform_arg);
    if (use_gall_precond && cg) {
        cg_opts_t prec_opts = main_cg_opts;
        prec_opts.precond_fn   = gall_precond_cb;
        prec_opts.precond_data = sg->gall;
        cg_state_set_opts(cg, &prec_opts);
    }

    tree_cov_matvec_inc_state_t inc;
    int inc_ok = (tree_cov_matvec_inc_state_init(&inc, sg->ts, n, sg->hap_to_indiv) == 0);

    if (!r || !r_tilde || !A || !B || !ws2_accum ||
        !gv || !cg || !inc_ok ||
        (use_fl && (!fl_log_ws2 || !fl_theta_in || !fl_theta_out))) {
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
        .gv = gv, .inc = &inc, .use_span = sg->uniform_arg,
    };
    hutch_opts_t hopts = {
        .nprobe = nprobe, .batch_size = 0, .use_gaussian = 0, .seed = 0,
        .probe_cg_tol = 0.0,        /* keep default 1e-3 */
        .probe_cg_max = cg_cap,     /* same bound as the main solve */
    };

    int uniform   = sg->uniform_arg;
    if (!silent)
        fprintf(stderr, "var_random: %s (n=%d, S=%d)\n",
                uniform ? "uniform-kinship REML"
                        : "per-tree REML EM",
                n, S);

    /* Per-tree ws2 bounds: floor prevents collapse to exact zero; ceiling
     * prevents runaway on a single tree.  In uniform mode the same ws2_max is
     * used as the cap on the single shared σ_u²/S value. */
    double ws2_max   = var_y / fmax(S, 1);
    double ws2_floor = 1e-12 * ws2_max;

    double sigma2_accum = 0.0;
    int    n_accum      = 0;
    double old_sigma2   = sg->sigma2;

    for (int iter = 0; iter < max_iter; iter++) {

        double A_0, B_0 = 0.0;
        vd.sigma2 = sg->sigma2;

        /* Rebuild G_all preconditioner from current ws2_mean + sigma2 (O(n)). */
        if (use_gall_precond) {
            double ws2_mean = 0.0;
            for (int s = 0; s < S; s++) ws2_mean += sg->ws2[s];
            ws2_mean /= fmax(S, 1);
            gall_precond_build(sg->gall, ws2_mean, sg->sigma2);
        }

        cg_solve(cg, v_matvec, &vd, r, r_tilde);

        A_0 = 0.0;
        for (int i = 0; i < n; i++) A_0 += r_tilde[i] * r_tilde[i];

        compute_quad_all(&inc, n, S, sg->S_orig, sg->prune_map,
                         r_tilde, A, gv, sg->uniform_arg);

        hutch_trace_all(cg, v_matvec, &vd, sg->ts, sg->hap_to_indiv,
                        n, S, sg->S_orig, sg->prune_map, B, &B_0, &hopts);

        /* REML EM update for segment weights.
         *
         * Per-tree mode (default): ws2[s] ← ws2[s] × A[s]/B[s].
         *   Fixed point A[s] = B[s] for every segment s.  Null segments shrink
         *   toward ws2_floor; active segments grow until their REML score is 0.
         *
         * Uniform mode: ws2[s] ≡ σ_u²/S for all s, where
         *   σ_u² ← σ_u² × Σ_s A[s] / Σ_s B[s].
         *   Equivalent to fitting a single random-effect variance against the
         *   uniformly-weighted kinship K = Σ_s G_s / S. */
        double ws2_mean = 0.0;
        if (uniform) {
            double A_sum = 0.0, B_sum = 0.0;
            for (int s = 0; s < S; s++) { A_sum += A[s]; B_sum += B[s]; }
            double ws2_old = sg->ws2[0];
            double ws2_new = ws2_old;
            if (B_sum > 1e-30) ws2_new = ws2_old * A_sum / B_sum;
            ws2_new = fmax(ws2_floor, fmin(ws2_max, ws2_new));
            for (int s = 0; s < S; s++) sg->ws2[s] = ws2_new;
            ws2_mean = ws2_new;
        } else {
            for (int s = 0; s < S; s++) {
                if (B[s] > 1e-30)
                    sg->ws2[s] *= A[s] / B[s];
                sg->ws2[s] = fmax(ws2_floor, fmin(ws2_max, sg->ws2[s]));
            }

            /* Block 2.6 — fused-lasso prox on log ws² (sparsity-on-levels +
             * sparsity-on-differences).  At lambda1=lambda2=0 (default) this
             * branch is skipped and behaviour is identical to the legacy EM. */
            if ((sg->fused_lasso_lambda1 > 0.0 || sg->fused_lasso_lambda2 > 0.0)
                && fl_log_ws2 != NULL) {
                double mu = 0.0;
                for (int s = 0; s < S; s++) {
                    fl_log_ws2[s] = log(sg->ws2[s]);
                    mu += fl_log_ws2[s];
                }
                mu /= S;
                for (int s = 0; s < S; s++)
                    fl_theta_in[s] = fl_log_ws2[s] - mu;

                fused_lasso_prox(fl_theta_in, fl_theta_out, S,
                                 sg->fused_lasso_lambda1,
                                 sg->fused_lasso_lambda2);

                for (int s = 0; s < S; s++) {
                    double w = exp(mu + fl_theta_out[s]);
                    sg->ws2[s] = fmax(ws2_floor, fmin(ws2_max, w));
                }
            }

            for (int s = 0; s < S; s++) ws2_mean += sg->ws2[s];
            ws2_mean /= S;
        }

        /* sigma2 REML EM update */
        if (B_0 > 0.0) {
            sg->sigma2 *= A_0 / B_0;
            sg->sigma2 = fmax(1e-10 * var_y, fmin(var_y, sg->sigma2));
        }

        double h2_rand = (var_y > 0.0) ?
            fmax(0.0, 1.0 - sg->sigma2 / var_y - sg->h2_snp) : 0.0;
        double rel_dsig = fabs(sg->sigma2 - old_sigma2) / (old_sigma2 + 1e-15);

        /* Polyak-Ruppert accumulation after burn-in */
        if (iter >= burn_in) {
            sigma2_accum += sg->sigma2;
            for (int s = 0; s < S; s++) ws2_accum[s] += sg->ws2[s];
            n_accum++;
        }

        if (!silent)
            fprintf(stderr, "var_random iter %3d: sigma2=%.4f  h2_rand=%.4f"
                    "  ws2_mean=%.3e  B0=%.2e  rel_dsig=%.2e%s\n",
                    iter + 1, sg->sigma2, h2_rand, ws2_mean, B_0, rel_dsig,
                    iter >= burn_in ? " [avg]" : "");

        old_sigma2 = sg->sigma2;
    }

    /* Apply Polyak-Ruppert averages */
    if (n_accum > 0) {
        sg->sigma2 = sigma2_accum / n_accum;
        for (int s = 0; s < S; s++)
            sg->ws2[s] = ws2_accum[s] / n_accum;
    }

    sg->h2_arg = (var_y > 0.0) ?
        fmax(0.0, 1.0 - sg->sigma2 / var_y - sg->h2_snp) : 0.0;
    sg->h2 = sg->h2_snp + sg->h2_arg;

    /* Post-hoc activation probability from converged ws2 values.
     *
     * After REML convergence, null trees cluster near ws2_floor and active
     * trees sit orders of magnitude higher.  We classify in log-space using a
     * threshold midway between floor and ceiling:
     *   ws2_thresh = sqrt(ws2_floor × ws2_max) = 1e-6 × ws2_max
     *
     * Score: log(ws2[s] / ws2_thresh), shifted by logit(phi_ref).
     *   phi_ref = phi_alpha / (phi_alpha + phi_beta)   (BB prior mean)
     *
     * Trees well above the threshold get rho ≈ 1; trees at the floor get
     * rho ≈ 0.  The Beta hyperprior (via phi_ref) shifts the decision
     * boundary: large phi_beta / phi_alpha (sparse prior, e.g. --pi0 0.9)
     * requires ws2[s] to be much larger than ws2_thresh to be called active. */
    {
        double ws2_thresh = sqrt(ws2_floor * ws2_max);   /* = 1e-6 × ws2_max */
        double phi_ref    = sg->phi_alpha / (sg->phi_alpha + sg->phi_beta + 1e-15);
        double logit_phi  = log(phi_ref / (1.0 - phi_ref + 1e-15));

        double rho_sum = 0.0;
        for (int s = 0; s < S; s++) {
            double lw     = log(fmax(sg->ws2[s], ws2_floor) / ws2_thresh);
            double score  = logit_phi + lw;
            score         = fmax(-SCORE_CAP, fmin(SCORE_CAP, score));
            double rho_s  = 1.0 / (1.0 + exp(-score));
            sg->rho[s * 2 + 0] = 1.0 - rho_s;  /* null */
            sg->rho[s * 2 + 1] = rho_s;          /* non-null */
            rho_sum += rho_s;
        }
        sg->phi_mix[sg->M - 1] = rho_sum / fmax(S, 1);
        /* Store the threshold scale so callers can recover tau2 ≈ ws2_thresh */
        sg->tau2_mix[sg->M - 1] = ws2_thresh;
    }

cleanup:
    cg_state_free(cg);
    if (inc_ok) tree_cov_matvec_inc_state_free(&inc);
    free(r);
    free(r_tilde);
    free(A);
    free(B);
    free(ws2_accum);
    free(gv);
    free(fl_log_ws2);
    free(fl_theta_in);
    free(fl_theta_out);
}

/* =========================================================================
 * var_random_emma: EMMA-based Stage 2 for uniform-ARG mode.
 *
 * Replaces var_random() when sg->uniform_arg=1 and sg->gall is set.
 * Instead of 100 EM iterations each costing ~700 tree-sequence passes, this
 * runs a 100-step golden-section search over the 1-D profile log-likelihood
 * in O(n²) time (one U^T y transform + 100×n evaluations).
 *
 * Falls back to var_random() when sg->gall is NULL.
 * ========================================================================= */
void var_random_emma(sargent_t *sg, int silent)
{
    if (!sg->gall) {
        if (!silent)
            fprintf(stderr, "var_random_emma: gall not available, "
                    "falling back to var_random\n");
        var_random(sg, silent);
        return;
    }

    int n = sg->n;
    int S = sg->S;

    /* Allocate ws2 if needed */
    if (!sg->ws2) {
        sg->ws2 = malloc((size_t)S * sizeof(double));
        if (!sg->ws2) { fprintf(stderr, "var_random_emma: OOM\n"); return; }
        double ws2_init = sg->tau2_mix[sg->M - 1] * 0.5;
        for (int s = 0; s < S; s++) sg->ws2[s] = ws2_init;
    }

    double var_y = vec_var(sg->y, n);

    if (!silent)
        fprintf(stderr, "var_random_emma: EMMA profile search (n=%d, S=%d)\n",
                n, S);

    /* Build residual r = y - X E[β] */
    double *r = malloc((size_t)n * sizeof(double));
    if (!r) { fprintf(stderr, "var_random_emma: OOM\n"); return; }
    memcpy(r, sg->y, (size_t)n * sizeof(double));
    if (sg->beta) {
        double *g = malloc((size_t)n * sizeof(double));
        if (!g) { free(r); fprintf(stderr, "var_random_emma: OOM\n"); return; }
        geno_rewind(sg->gr);
        for (int j = 0; j < sg->p; j++) {
            if (geno_read_variant(sg->gr, g, GENO_STANDARDIZE) <= 0) break;
            if (sg->beta[j] == 0.0) continue;
            for (int i = 0; i < n; i++) r[i] -= sg->beta[j] * g[i];
        }
        free(g);
    }

    /* EMMA profile-likelihood search for (ws2_shared, sigma2) */
    double var_r = vec_var(r, n);
    double ws2_shared, sigma2;
    int niters = gall_emma_reml(sg->gall, r, &ws2_shared, &sigma2);
    free(r);

    /* Boundary check: if sigma2 << var(r), EMMA hit the flat ML boundary
     * (G_all is full-rank and absorbed all residual variance).  This happens
     * at small n or near-zero h²_arg.  Fall back to regularised EM. */
    if (sigma2 < 1e-3 * var_r) {
        if (!silent)
            fprintf(stderr, "var_random_emma: boundary (sigma2=%.2e << var_r=%.2e),"
                    " falling back to var_random\n", sigma2, var_r);
        var_random(sg, silent);
        return;
    }

    /* Clamp to physically sensible range */
    double ws2_max   = var_y / fmax(S, 1);
    double ws2_floor = 1e-12 * ws2_max;
    ws2_shared = fmax(ws2_floor, fmin(ws2_max, ws2_shared));
    sigma2     = fmax(1e-10 * var_y, fmin(var_y, sigma2));

    /* Apply: uniform ws2 across all segments */
    for (int s = 0; s < S; s++) sg->ws2[s] = ws2_shared;
    sg->sigma2 = sigma2;

    sg->h2_arg = (var_y > 0.0) ?
        fmax(0.0, 1.0 - sg->sigma2 / var_y - sg->h2_snp) : 0.0;
    sg->h2 = sg->h2_snp + sg->h2_arg;

    double h2_rand = fmax(0.0, 1.0 - sigma2 / var_y - sg->h2_snp);
    if (!silent)
        fprintf(stderr, "var_random_emma: done in %d steps  "
                "ws2=%.4g  sigma2=%.4f  h2_rand=%.4f\n",
                niters, ws2_shared, sigma2, h2_rand);

    /* Post-hoc activation rho: uniform mode → all segments equal, so use a
     * single global ws2_thresh derived from the floor/ceiling. */
    {
        double ws2_thresh = sqrt(ws2_floor * ws2_max);
        double phi_ref    = sg->phi_alpha / (sg->phi_alpha + sg->phi_beta + 1e-15);
        double logit_phi  = log(phi_ref / (1.0 - phi_ref + 1e-15));

        double lw    = log(fmax(ws2_shared, ws2_floor) / ws2_thresh);
        double score = fmax(-SCORE_CAP, fmin(SCORE_CAP, logit_phi + lw));
        double rho_s = 1.0 / (1.0 + exp(-score));

        for (int s = 0; s < S; s++) {
            sg->rho[s * 2 + 0] = 1.0 - rho_s;
            sg->rho[s * 2 + 1] = rho_s;
        }
        sg->phi_mix[sg->M - 1] = rho_s;
        sg->tau2_mix[sg->M - 1] = ws2_thresh;
    }
}

/* =========================================================================
 * var_joint: interleaved coordinate descent for β, ws², σ² (joint EM)
 * =========================================================================
 *
 * Each outer iteration performs:
 *   1. r = y - Xβ           (residual)
 *   2. r̃ = V^{-1} r         (one CG solve; V = Σ ws²·G_s + σ²I)
 *   3. Compute A[s], B[s], A_0, B_0   (per-segment quadratic forms + traces)
 *   4. ONE sweep over β_j using V-aware precision v_jj ≈ tr(V^{-1}) = B_0
 *      and t_j = x_j^T r̃ + β_j_old · v_jj  (leave-one-out add-back)
 *   5. Update π_k (BayesR mixture proportions)
 *   6. Update ws²[s] via REML EM (per-tree or uniform mode)
 *   7. Update σ² via REML EM (multiplicative)
 *
 * The β update uses the current (stale w.r.t. β-sweep) r̃; r̃ is recomputed at
 * the start of the next iter.  This is the standard mean-field VB pattern:
 * each parameter sees the latest value of the others between outer iterations.
 *
 * Initialised with β = 0 and warm uniform ws² so the very first β-sweep
 * already sees a non-trivial random offset and cannot cannibalise.
 *
 * Post-hoc rho is computed after convergence, matching var_random.
 * ====================================================================== */
void var_joint(sargent_t *sg, int silent)
{
    int n      = sg->n;
    int p      = sg->p;
    int S      = sg->S;
    int K      = sg->K;
    int max_outer = sg->em_iters > 0 ? sg->em_iters : 100;

    /* Allocate β (zero) if not already done */
    if (!sg->beta) {
        sg->beta = calloc((size_t)p, sizeof(double));
        if (!sg->beta) { fprintf(stderr, "var_joint: out of memory\n"); return; }
    }

    /* Allocate ws² with warm uniform init if not already done */
    if (!sg->ws2) {
        sg->ws2 = malloc((size_t)S * sizeof(double));
        if (!sg->ws2) { fprintf(stderr, "var_joint: out of memory\n"); return; }
        double ws2_init = sg->tau2_mix[sg->M - 1] * 0.5;
        for (int s = 0; s < S; s++) sg->ws2[s] = ws2_init;
    }

    /* Working buffers */
    double   *r          = malloc((size_t)n * sizeof(double));
    double   *r_tilde    = calloc((size_t)n,  sizeof(double));
    double   *A          = malloc((size_t)S * sizeof(double));
    double   *B          = malloc((size_t)S * sizeof(double));
    double   *g          = malloc((size_t)n * sizeof(double));
    double   *log_g      = malloc((size_t)K * sizeof(double));
    double   *gamma      = malloc((size_t)K * sizeof(double));
    double   *gamma_sum  = malloc((size_t)K * sizeof(double));

    double    *gv        = malloc((size_t)n * sizeof(double));

    cg_state_t *cg = cg_state_new(n, NULL);

    tree_cov_matvec_inc_state_t inc;
    int inc_ok = (tree_cov_matvec_inc_state_init(&inc, sg->ts, n, sg->hap_to_indiv) == 0);

    if (!r || !r_tilde || !A || !B || !g || !log_g || !gamma ||
        !gamma_sum || !gv || !cg || !inc_ok) {
        fprintf(stderr, "var_joint: out of memory\n");
        goto cleanup;
    }

    double var_y     = vec_var(sg->y, n);
    double ws2_max   = var_y / fmax(S, 1);
    double ws2_floor = 1e-12 * ws2_max;

    /* Phase 4.5 divergence-detection snapshots: at the start of every iter
     * we snapshot (β, ws², σ²); if the iter produced a wild rel_dsig we
     * restore to break out of an exploding polish trajectory. */
    double *beta_snap   = malloc((size_t)p * sizeof(double));
    double *ws2_snap    = malloc((size_t)S * sizeof(double));
    double  sigma2_snap = 0.0;
    double  h2_snp_snap = 0.0;
    double  h2_arg_snap = 0.0;
    if (!beta_snap || !ws2_snap) {
        fprintf(stderr, "var_joint: out of memory (snapshot)\n");
        free(beta_snap); free(ws2_snap);
        goto cleanup;
    }

    v_matvec_data_t vd = {
        .ts = sg->ts, .ws2 = sg->ws2, .sigma2 = sg->sigma2,
        .n = n, .S = S, .S_orig = sg->S_orig, .prune_map = sg->prune_map,
        .hap_to_indiv = sg->hap_to_indiv,
        .gv = gv, .inc = &inc, .use_span = sg->uniform_arg,
    };
    hutch_opts_t hopts = {
        .nprobe = sg->nprobe, .batch_size = 0, .use_gaussian = 0, .seed = 0,
    };

    int uniform   = sg->uniform_arg;
    if (!silent)
        fprintf(stderr,
                "var_joint: interleaved coordinate descent "
                "(n=%d, p=%d, S=%d, %s, max_outer=%d)\n",
                n, p, S,
                uniform ? "uniform-random" : "per-tree",
                max_outer);

    double old_sigma2   = sg->sigma2;

    for (int outer = 0; outer < max_outer; outer++) {

        /* --- 0. Snapshot (β, ws², σ², h²_snp, h²_arg) for divergence-revert --- */
        memcpy(beta_snap, sg->beta, (size_t)p * sizeof(double));
        memcpy(ws2_snap,  sg->ws2,  (size_t)S * sizeof(double));
        sigma2_snap = sg->sigma2;
        h2_snp_snap = sg->h2_snp;
        h2_arg_snap = sg->h2_arg;

        /* --- 1. Compute residual r = y - X β --- */
        memcpy(r, sg->y, (size_t)n * sizeof(double));
        int any_beta = 0;
        for (int j = 0; j < p; j++) {
            if (sg->beta[j] != 0.0) { any_beta = 1; break; }
        }
        if (any_beta) {
            geno_rewind(sg->gr);
            for (int j = 0; j < p; j++) {
                if (geno_read_variant(sg->gr, g, GENO_STANDARDIZE) <= 0) break;
                if (sg->beta[j] == 0.0) continue;
                for (int i = 0; i < n; i++) r[i] -= sg->beta[j] * g[i];
            }
        }

        /* --- 2. r̃ = V^{-1} r --- */
        vd.sigma2 = sg->sigma2;
        cg_solve(cg, v_matvec, &vd, r, r_tilde);

        /* --- 3. A_0, A[s], B_0, B[s] --- */
        double A_0 = 0.0, B_0 = 0.0;
        for (int i = 0; i < n; i++) A_0 += r_tilde[i] * r_tilde[i];
        compute_quad_all(&inc, n, S, sg->S_orig, sg->prune_map,
                         r_tilde, A, gv, sg->uniform_arg);
        hutch_trace_all(cg, v_matvec, &vd, sg->ts, sg->hap_to_indiv,
                        n, S, sg->S_orig, sg->prune_map, B, &B_0, &hopts);

        /* --- 4. One β sweep with constant V-aware precision v_jj ≈ tr(V^{-1}) ---
         *
         * The exact variational update would use v_jj = x_j^T V^{-1} x_j, but
         * computing that per SNP requires p CG solves per outer iter — not
         * feasible at biobank scale (p ~ 1M).  We use the constant average
         * v_jj = tr(V^{-1}) ≈ B_0 for all SNPs.  This matches BOLT-LMM's
         * approximation (see Loh et al. 2015): rely on iteration with the
         * V-aware residual r̃ to capture the V-correction globally, not
         * per-SNP.  Known trade-off in our setting: when X and K share
         * substantial signal (as for ARG-based kinship), this approximation
         * over-shrinks β and pushes variance into the random component.
         * Joint EM is therefore a good predictor but a poor variance-component
         * estimator — use sequential SARGENT for the latter. */
        double v_jj = (B_0 > 0.0) ? B_0 : (double)n / sg->sigma2;
        memset(gamma_sum, 0, (size_t)K * sizeof(double));
        double max_dbeta = 0.0;
        int    n_read    = 0;

        geno_rewind(sg->gr);
        for (int j = 0; j < p; j++) {
            if (geno_read_variant(sg->gr, g, GENO_STANDARDIZE) <= 0) break;
            n_read++;
            double beta_old = sg->beta[j];

            /* t_j = x_j^T r̃ + β_old · v_jj  (LOO add-back via v_jj approximation) */
            double xrt = 0.0;
            for (int i = 0; i < n; i++) xrt += g[i] * r_tilde[i];
            double t_j = xrt + beta_old * v_jj;

            /* BayesR mixture posterior with V-aware precision v_jj */
            double z_sq = t_j * t_j / v_jj;
            double lmax = -1e300;
            for (int k = 0; k < K; k++) {
                double sk2 = sg->sigma2_mix[k];
                double log_pi_k = log(sg->pi_mix[k] + PI_FLOOR);
                if (sk2 == 0.0) {
                    log_g[k] = log_pi_k - 0.5 * z_sq;
                } else {
                    double s_k = v_jj * sk2 + 1.0;
                    log_g[k] = log_pi_k - 0.5 * (log(s_k) + z_sq / s_k);
                }
                if (log_g[k] > lmax) lmax = log_g[k];
            }
            double sum_g = 0.0;
            for (int k = 0; k < K; k++) {
                gamma[k] = exp(log_g[k] - lmax);
                sum_g += gamma[k];
            }
            for (int k = 0; k < K; k++) gamma[k] /= sum_g;

            double beta_raw = 0.0;
            for (int k = 0; k < K; k++) {
                if (gamma[k] == 0.0 || sg->sigma2_mix[k] == 0.0) continue;
                double v_post = 1.0 / (v_jj + 1.0 / sg->sigma2_mix[k]);
                beta_raw += gamma[k] * t_j * v_post;
            }
            /* Damped β update: β_new = (1-ω)·β_old + ω·β_raw, ω=0.5.
             *
             * Phase 4.5: r̃ is held stale through the entire β-sweep (line 829),
             * so with p ≫ n_train and SNPs in LD all SNPs see the same residual
             * and cumulatively overshoot.  Halving the step kills the
             * cumulative drift without changing the EM fixed point. */
            const double EM_DAMP_BETA = 0.5;
            double beta_new = (1.0 - EM_DAMP_BETA) * beta_old + EM_DAMP_BETA * beta_raw;
            sg->beta[j] = beta_new;
            for (int k = 0; k < K; k++) gamma_sum[k] += gamma[k];

            double delta = beta_new - beta_old;
            if (fabs(delta) > max_dbeta) max_dbeta = fabs(delta);

            /* Update r so subsequent β_j sweeps see the change.
             * r̃ is NOT updated within the sweep — refreshed at start of next iter. */
            for (int i = 0; i < n; i++) r[i] -= delta * g[i];
        }

        /* --- 5. Update mixture proportions --- */
        if (n_read > 0)
            for (int k = 0; k < K; k++)
                sg->pi_mix[k] = gamma_sum[k] / n_read;

        /* --- 6. ws²[s] update (per-tree REML EM or uniform mode) --- */
        double ws2_mean = 0.0;
        if (uniform) {
            double A_sum = 0.0, B_sum = 0.0;
            for (int s = 0; s < S; s++) { A_sum += A[s]; B_sum += B[s]; }
            double ws2_old = sg->ws2[0];
            double ws2_new = ws2_old;
            if (B_sum > 1e-30) ws2_new = ws2_old * A_sum / B_sum;
            ws2_new = fmax(ws2_floor, fmin(ws2_max, ws2_new));
            for (int s = 0; s < S; s++) sg->ws2[s] = ws2_new;
            ws2_mean = ws2_new;
        } else {
            for (int s = 0; s < S; s++) {
                if (B[s] > 1e-30) sg->ws2[s] *= A[s] / B[s];
                sg->ws2[s] = fmax(ws2_floor, fmin(ws2_max, sg->ws2[s]));
                ws2_mean += sg->ws2[s];
            }
            ws2_mean /= S;
        }

        /* --- 7. σ² update (Phase 4.5: damped + ratio-clipped + raised floor) ---
         *
         * The plain multiplicative update σ² *= A_0/B_0 oscillates wildly
         * when V is near-singular (which it becomes as σ² shrinks and ws²
         * is small).  Three layers of stabilisation:
         *
         *   1. Damp with ω=0.25: σ²_new = σ²·(1 + ω·(A_0/B_0 - 1)).
         *   2. Clip per-iter step ratio to [0.5, 2.0].  Even with
         *      damping, A_0/B_0 can swing from 0.3 to 3+ in a single iter
         *      when CG produces near-singular r_tilde; the cap prevents
         *      single-step blow-ups.
         *   3. Floor σ² at 0.001·var_y (was 1e-10·var_y).  Below 0.001·var_y,
         *      V is dominated by the random component (ws² >> σ²) and CG
         *      is ill-conditioned, which is exactly what drives the
         *      A_0/B_0 swings.  Keeping σ² above the floor keeps V away
         *      from singularity.
         *
         * Note: differs from var_fixed's asymmetric (full-step descent /
         * half-step ascent) σ² damping.  Both target the same root cause
         * (BayesR mixture-snap cascade) but var_joint's Newton-ratio
         * proposal is unbounded above so it needs the symmetric+clip form.
         * See ROADMAP "σ² damping convention" subsection. */
        const double EM_DAMP_SIGMA = 0.25;
        if (B_0 > 0.0) {
            double ratio = A_0 / B_0;
            double step  = 1.0 + EM_DAMP_SIGMA * (ratio - 1.0);
            step         = fmax(0.5, fmin(2.0, step));
            sg->sigma2  *= step;
            sg->sigma2   = fmax(1e-3 * var_y, fmin(var_y, sg->sigma2));
        }

        /* --- 8. h² accounting (Phase 4.5: report var(E[u])/var(y) properly) ---
         *
         * BLUP of u given y is E[u | y, β] = K_g V^{-1} r = (V-σ²I)V^{-1}r
         *                                  = r - σ²·r_tilde.
         * h²_snp = var(Xβ)/var(y) = (var_y - var_r)/var_y    (since y centered)
         * h²_arg = var(E[u])/var(y)
         * h²     = h²_snp + h²_arg
         *
         * This replaces the chunk-4 reporting h²_arg = (var_r - σ²)/var_y,
         * which treated h²_arg as "everything not in X β minus noise" --
         * an over-estimate when E[u] doesn't actually fit (the Nurdan
         * σ²-collapse failure mode). */
        double var_r = vec_var(r, n);
        double var_eu = 0.0;
        {
            double mean_eu = 0.0;
            for (int i = 0; i < n; i++) mean_eu += r[i] - sg->sigma2 * r_tilde[i];
            mean_eu /= n;
            for (int i = 0; i < n; i++) {
                double eu = (r[i] - sg->sigma2 * r_tilde[i]) - mean_eu;
                var_eu += eu * eu;
            }
            var_eu /= n;
        }
        sg->h2_snp = (var_y > 0.0) ? fmax(0.0, 1.0 - var_r / var_y) : 0.0;
        sg->h2_arg = (var_y > 0.0) ? fmax(0.0, fmin(1.0, var_eu / var_y)) : 0.0;
        sg->h2     = sg->h2_snp + sg->h2_arg;

        double rel_dsig = fabs(sg->sigma2 - old_sigma2) / (old_sigma2 + 1e-15);

        if (!silent)
            fprintf(stderr,
                    "var_joint iter %3d: max_dβ=%.2e  sigma2=%.4f  h2_f=%.4f  h2_r=%.4f"
                    "  ws2_mean=%.3e  v_jj=%.2e  rel_dsig=%.2e\n",
                    outer + 1, max_dbeta, sg->sigma2, sg->h2_snp, sg->h2_arg,
                    ws2_mean, v_jj, rel_dsig);

        /* --- 9. Phase 4.5: divergence detection --------------------------
         *
         * The polish trajectory should track a smooth EM descent.  rel_dsig
         * above 0.05 (5% per-iter relative σ² change) signals the start of
         * an oscillation; observed empirically on the s42_balanced debug
         * case where iter 9 jumped to rel_dsig=0.21 and the polish then
         * diverged.  Revert to the start-of-iter snapshot and stop.
         *
         * Skip the first iter (cold-start can produce a large legitimate
         * step from the warm-start init to the polish minimum). */
        if (outer >= 1 && rel_dsig > 0.05) {
            if (!silent)
                fprintf(stderr,
                        "var_joint: divergence detected at iter %d "
                        "(rel_dsig=%.2e > 0.05); reverting to iter %d state\n",
                        outer + 1, rel_dsig, outer);
            memcpy(sg->beta, beta_snap, (size_t)p * sizeof(double));
            memcpy(sg->ws2,  ws2_snap,  (size_t)S * sizeof(double));
            sg->sigma2 = sigma2_snap;
            sg->h2_snp = h2_snp_snap;
            sg->h2_arg = h2_arg_snap;
            sg->h2     = h2_snp_snap + h2_arg_snap;
            break;
        }

        old_sigma2 = sg->sigma2;
    }

    /* Final h2 values reflect the last iteration's coherent (β, ws², σ²)
     * triple; var_r came from the very last β sweep.  All quantities written
     * to .beta.tsv and .w2.tsv (and these h² fields) are from the same
     * iteration — no PR averaging — so downstream consumers can trust
     * h²_fixed + h²_random + σ²/var_y ≈ 1 without reconstruction. */

    /* --- Post-hoc rho (identical to var_random) --- */
    {
        double ws2_thresh = sqrt(ws2_floor * ws2_max);
        double phi_ref    = sg->phi_alpha / (sg->phi_alpha + sg->phi_beta + 1e-15);
        double logit_phi  = log(phi_ref / (1.0 - phi_ref + 1e-15));
        double rho_sum    = 0.0;
        for (int s = 0; s < S; s++) {
            double lw    = log(fmax(sg->ws2[s], ws2_floor) / ws2_thresh);
            double score = logit_phi + lw;
            score        = fmax(-SCORE_CAP, fmin(SCORE_CAP, score));
            double rho_s = 1.0 / (1.0 + exp(-score));
            sg->rho[s * 2 + 0] = 1.0 - rho_s;
            sg->rho[s * 2 + 1] = rho_s;
            rho_sum += rho_s;
        }
        sg->phi_mix[sg->M - 1] = rho_sum / fmax(S, 1);
        sg->tau2_mix[sg->M - 1] = ws2_thresh;
    }

cleanup:
    cg_state_free(cg);
    if (inc_ok) tree_cov_matvec_inc_state_free(&inc);
    free(beta_snap); free(ws2_snap);
    free(r); free(r_tilde); free(A); free(B); free(g);
    free(log_g); free(gamma); free(gamma_sum);
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
