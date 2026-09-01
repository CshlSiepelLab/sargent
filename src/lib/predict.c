/*
 * SARGENT: Sparse ARG-based ENcoding of Traits
 *
 * Copyright (c) 2026, Cold Spring Harbor Laboratory
 * All rights reserved.
 *
 * This file is part of SARGENT and is distributed under the BSD 3-Clause License.
 * See the LICENSE file in the project root for details.
 *
 * predict.c — native out-of-sample prediction (full five-arm ladder).
 * Reproduces analysis/predict.py; see test/eval/native_predict_design.md for
 * the exact math and the < 1e-6 merge gate.
 */

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <time.h>

/* ---- per-arm wall clock ---------------------------------------------------
 * Every invocation computes EVERY arm -- --snp-only, --uniform-arg and --joint
 * all emit snp_lmm, sargent_s, sargent_a and sargent -- so the process wall clock
 * cannot be attributed to any one flavour, and differencing whole runs is both
 * imprecise and easy to misread (it was: 186 s was quoted as SARGENT-UR's cost
 * when the ARG-LMM fit is 0.2 s of it).  These counters time each arm's own
 * fit, so the JSON carries an honest per-flavour cost.  Emit-only: no arm's
 * numbers change. */
/* Phase timers filled in by the driver (sargent.c) for work that happens
 * BEFORE predict_run is called, so the JSON can attribute it. */
double sg_t_stage1 = 0.0, sg_t_stage1_f = 0.0, sg_t_gall = 0.0,
       sg_t_joint = 0.0, sg_t_stage2 = 0.0;

double arm_now(void)
{
    struct timespec ts;
    clock_gettime(CLOCK_MONOTONIC, &ts);
    return (double)ts.tv_sec + 1e-9 * (double)ts.tv_nsec;
}

#include <tskit.h>
#include "predict.h"
#include "geno.h"
#include "sargent.h"
#include "branchsel.h"

/* Which arms this invocation computes.  SG_MODEL_ALL is the default and
 * reproduces the historical behaviour exactly. */
static unsigned g_arm_mask = SG_MODEL_ALL;
void predict_set_models(unsigned mask) { g_arm_mask = mask ? mask : SG_MODEL_ALL; }
#define ARM_ON(bits) ((g_arm_mask & (unsigned)(bits)) != 0u)

/* ------------------------------------------------------------------------- *
 * Portable LAPACK declarations (Fortran ABI; same as gall.c's dsyev_).
 * dpotrf: Cholesky factor; dpotri: inverse from factor; dposv: SPD solve.
 * ------------------------------------------------------------------------- */
extern void dpotrf_(const char *uplo, const int *n, double *a, const int *lda,
                    int *info);
extern void dpotri_(const char *uplo, const int *n, double *a, const int *lda,
                    int *info);
extern void dposv_(const char *uplo, const int *n, const int *nrhs,
                   double *a, const int *lda, double *b, const int *ldb,
                   int *info);
extern void dsyev_(const char *jobz, const char *uplo, const int *n, double *a,
                   const int *lda, double *w, double *work, const int *lwork,
                   int *info);


/* Realized-variance h^2 for a single-kernel BLUP (SARGENT-UR / sargent_a).
 *
 * The joint chain reports h^2 as Var(Xb + u) / [Var(Xb + u) + sigma2_e] (BSLMM
 * eq. 13) and is nearly unbiased for it; the REML plug-in su2/(su2+se2) used for
 * this arm climbs with phi (+0.028 -> +0.094 at h2=0.8, +0.051 -> +0.123 at
 * h2=0.5) because su2 inflates to absorb a misspecified kernel.  Leaving the two
 * SARGENT arms on different estimators put two lines in one panel that behave
 * differently for no reason a reader could infer.
 *
 * Same estimator, computed analytically since this arm has no sampler: in the
 * eigenbasis of the centred K, with r = U^T (y - mean),
 *     m_i = su2*lam_i/(su2*lam_i + se2) * r_i        posterior mean
 *     v_i = su2*lam_i*se2/(su2*lam_i + se2)          posterior variance
 * and Var(u) = (1/n) sum (m_i^2 + v_i) -- the mean-plus-variance form, not the
 * variance of the posterior mean, which would under-count by exactly the
 * shrinkage.  One dsyev at n = 1000 costs ~1 s against this arm's ~188 s.
 * Returns -1 if the decomposition fails, leaving the caller's value untouched. */
static int realized_h2_1k(const double *K, const double *y, int n,
                          double su2, double se2, double *out_h2)
{
    double *A = malloc((size_t)n * n * sizeof(double));
    double *w = malloc((size_t)n * sizeof(double));
    double *r = malloc((size_t)n * sizeof(double));
    if (!A || !w || !r) { free(A); free(w); free(r); return -1; }
    memcpy(A, K, (size_t)n * n * sizeof(double));

    int info = 0, lwork = -1; double qwork = 0.0;
    dsyev_("V", "U", &n, A, &n, w, &qwork, &lwork, &info);
    if (info != 0) { free(A); free(w); free(r); return -1; }
    lwork = (int)qwork;
    double *work = malloc((size_t)lwork * sizeof(double));
    if (!work) { free(A); free(w); free(r); return -1; }
    dsyev_("V", "U", &n, A, &n, w, work, &lwork, &info);
    free(work);
    if (info != 0) { free(A); free(w); free(r); return -1; }

    double mean = 0.0;
    for (int i = 0; i < n; i++) mean += y[i];
    mean /= (double)n;
    for (int k = 0; k < n; k++) {          /* r = U^T (y - mean) */
        double s = 0.0;
        for (int i = 0; i < n; i++) s += A[(size_t)k * n + i] * (y[i] - mean);
        r[k] = s;
    }
    double vg = 0.0;
    for (int k = 0; k < n; k++) {
        double lam = w[k] > 0.0 ? w[k] : 0.0;
        double d = su2 * lam + se2; if (d < 1e-12) d = 1e-12;
        double m = (su2 * lam / d) * r[k];
        double v = (su2 * lam) * se2 / d;
        vg += m * m + v;
    }
    vg /= (double)n;
    free(A); free(w); free(r);
    if (vg + se2 <= 1e-12) return -1;
    *out_h2 = vg / (vg + se2);
    return 0;
}

/* ------------------------------------------------------------------------- *
 * Small helpers.
 * ------------------------------------------------------------------------- */

/* Symmetrize a full n×n matrix in place from its lower triangle (col-major),
 * which is the upper triangle in row-major.  We store everything row-major and
 * treat symmetric matrices as such, so this mirrors triangle a[i][j]<->a[j][i]. */
static void symmetrize_full(double *A, int n)
{
    for (int i = 0; i < n; i++)
        for (int j = i + 1; j < n; j++) {
            double v = A[(size_t)i * n + j];
            A[(size_t)j * n + i] = v;
        }
}

/* out[i] = (A @ x)[i] for row-major n×n symmetric A and length-n x. */
static void matvec(const double *A, const double *x, double *out, int n)
{
    for (int i = 0; i < n; i++) {
        const double *Ai = A + (size_t)i * n;
        double s = 0.0;
        for (int j = 0; j < n; j++) s += Ai[j] * x[j];
        out[i] = s;
    }
}

/* JSON number formatting matching Python json.dump (NaN/Infinity literals). */
static void json_num(FILE *f, double v)
{
    if (isnan(v))            fputs("NaN", f);
    else if (isinf(v))       fputs(v > 0 ? "Infinity" : "-Infinity", f);
    else                     fprintf(f, "%.12g", v);
}

/* ------------------------------------------------------------------------- *
 * Metrics (predict.py:metrics).
 * ------------------------------------------------------------------------- */
typedef struct { double r, r2, r2_calib, mse; } metrics_t;

static metrics_t compute_metrics(const double *yhat, const double *y, int n)
{
    metrics_t m;
    double my = 0.0;
    for (int i = 0; i < n; i++) my += y[i];
    my /= n;
    double ss_res = 0.0, ss_tot = 0.0;
    for (int i = 0; i < n; i++) {
        double d = y[i] - yhat[i];
        ss_res += d * d;
        double t = y[i] - my;
        ss_tot += t * t;
    }
    m.r2 = (ss_tot > 0.0) ? 1.0 - ss_res / ss_tot : NAN;
    m.mse = ss_res / n;

    /* Pearson r */
    double mh = 0.0;
    for (int i = 0; i < n; i++) mh += yhat[i];
    mh /= n;
    double sxy = 0.0, sxx = 0.0, syy = 0.0;
    for (int i = 0; i < n; i++) {
        double a = yhat[i] - mh, b = y[i] - my;
        sxy += a * b; sxx += a * a; syy += b * b;
    }
    if (sxx > 0.0 && syy > 0.0)
        m.r = sxy / sqrt(sxx * syy);
    else
        m.r = NAN;
    m.r2_calib = (m.r == m.r) ? m.r * m.r : NAN;
    return m;
}

/* ------------------------------------------------------------------------- *
 * center_kinship (predict.py:center_kinship).
 * Centers K_train (H K H) and, consistently, the cross block; then mean-diag=1.
 * Both matrices are modified in place.  K_cross may be NULL.
 * ------------------------------------------------------------------------- */
static void center_kinship(double *K_train, int n_tr,
                           double *K_cross, int n_cr)
{
    double *row = malloc((size_t)n_tr * sizeof(double)); /* rowmean = colmean (sym) */
    double grand = 0.0;
    for (int i = 0; i < n_tr; i++) {
        const double *Ki = K_train + (size_t)i * n_tr;
        double s = 0.0;
        for (int j = 0; j < n_tr; j++) s += Ki[j];
        row[i] = s / n_tr;
        grand += s;
    }
    grand /= (double)n_tr * n_tr;

    for (int i = 0; i < n_tr; i++) {
        double *Ki = K_train + (size_t)i * n_tr;
        for (int j = 0; j < n_tr; j++)
            Ki[j] = Ki[j] - row[i] - row[j] + grand;
    }
    double d = 0.0;
    for (int i = 0; i < n_tr; i++) d += K_train[(size_t)i * n_tr + i];
    d /= n_tr;

    if (d > 0.0) {
        for (size_t k = 0; k < (size_t)n_tr * n_tr; k++) K_train[k] /= d;
    }

    if (K_cross) {
        /* cross uses TRAIN col means (== row[]) and grand; its own row means. */
        for (int c = 0; c < n_cr; c++) {
            double *Kc = K_cross + (size_t)c * n_tr;
            double s = 0.0;
            for (int j = 0; j < n_tr; j++) s += Kc[j];
            double rmean = s / n_tr;
            for (int j = 0; j < n_tr; j++)
                Kc[j] = Kc[j] - rmean - row[j] + grand;
            if (d > 0.0)
                for (int j = 0; j < n_tr; j++) Kc[j] /= d;
        }
    }
    free(row);
}

/* ------------------------------------------------------------------------- *
 * fit_blup_m (predict.py:fit_blup, generalized to nK variance components).
 *
 * V = sum_k su2[k] K[k] + s2 I, fitted by the same multiplicative REML-EM as
 * the single-component original, then
 *
 *     alpha = A^{-1} y_c,  A = sum_k (su2[k]/su2[0]) K[k] + (s2/su2[0]) I
 *
 * so that yhat = sum_k (su2[k]/su2[0]) K_cross[k] alpha.  That scaling is
 * chosen so nK = 1 collapses to the original expression (K + lam I)^{-1} y_c
 * with lam = s2/su2 term for term, keeping every pre-existing arm bit-identical.
 *
 * K[k] are row-major n x n (centered).  P_out, when non-NULL, receives the REML
 * projector P = V^{-1} - V^{-1} 1 1' V^{-1} / (1' V^{-1} 1) at the fitted
 * variance components — the background against which branch scores are taken.
 * ratio_out, when non-NULL, receives su2[k]/su2[0] for the cross-block combine.
 * ------------------------------------------------------------------------- */
static int fit_blup_m(const double *y, const double *const *K, int nK, int n,
                      double *out_s2, double *out_su2,
                      double *alpha, double *out_mean_y,
                      double *ratio_out, double *P_out)
{
    const int max_iter = 50;
    const double tol = 1e-5;

    double mean = 0.0;
    for (int i = 0; i < n; i++) mean += y[i];
    mean /= n;
    double *yc = malloc((size_t)n * sizeof(double));
    double var_y = 0.0;
    for (int i = 0; i < n; i++) { yc[i] = y[i] - mean; var_y += yc[i] * yc[i]; }
    var_y /= n;

    double su2[8], s2 = 0.5 * var_y;
    if (nK < 1 || nK > 8) { free(yc); return -1; }
    for (int k = 0; k < nK; k++) su2[k] = 0.5 * var_y / nK;

    double *V     = malloc((size_t)n * n * sizeof(double));
    double *Vinv  = malloc((size_t)n * n * sizeof(double));
    double *Vy    = malloc((size_t)n * sizeof(double));
    double *KVy   = malloc((size_t)n * sizeof(double));
    if (!yc || !V || !Vinv || !Vy || !KVy) {
        free(yc); free(V); free(Vinv); free(Vy); free(KVy);
        return -1;
    }

    for (int iter = 0; iter < max_iter; iter++) {
        /* V = sum_k su2[k]*K[k] + s2*I */
        for (size_t q = 0; q < (size_t)n * n; q++) V[q] = su2[0] * K[0][q];
        for (int k = 1; k < nK; k++)
            for (size_t q = 0; q < (size_t)n * n; q++) V[q] += su2[k] * K[k][q];
        for (int i = 0; i < n; i++) V[(size_t)i * n + i] += s2;

        /* Cholesky; on failure add 1e-6 I and retry once (matches predict.py). */
        memcpy(Vinv, V, (size_t)n * n * sizeof(double));
        int info = 0;
        dpotrf_("L", &n, Vinv, &n, &info);
        if (info != 0) {
            for (int i = 0; i < n; i++) V[(size_t)i * n + i] += 1e-6;
            memcpy(Vinv, V, (size_t)n * n * sizeof(double));
            dpotrf_("L", &n, Vinv, &n, &info);
            if (info != 0) { free(yc); free(V); free(Vinv); free(Vy); free(KVy); return -1; }
        }
        /* Vinv = V^{-1} from the Cholesky factor. */
        dpotri_("L", &n, Vinv, &n, &info);
        if (info != 0) { free(yc); free(V); free(Vinv); free(Vy); free(KVy); return -1; }
        symmetrize_full(Vinv, n);

        matvec(Vinv, yc, Vy, n);           /* Vinv_y */

        double a_e = 0.0;
        for (int i = 0; i < n; i++) a_e += Vy[i] * Vy[i];
        double tr_Vinv = 0.0;
        for (int i = 0; i < n; i++) tr_Vinv += Vinv[(size_t)i * n + i];

        double floor_v = 1e-8 * var_y;
        int converged = 1;
        for (int k = 0; k < nK; k++) {
            matvec(K[k], Vy, KVy, n);      /* K_k Vinv_y */
            double a_u = 0.0;
            for (int i = 0; i < n; i++) a_u += KVy[i] * Vy[i];

            double tr_VinvK = 0.0;
            for (int i = 0; i < n; i++) {
                const double *Vi = Vinv + (size_t)i * n;
                const double *Ki = K[k] + (size_t)i * n;
                double row = 0.0;
                for (int j = 0; j < n; j++) row += Vi[j] * Ki[j];
                tr_VinvK += row;
            }

            double new_su2 = su2[k] * a_u / fmax(tr_VinvK, 1e-30);
            if (new_su2 < floor_v) new_su2 = floor_v;
            if (!(fabs(new_su2 - su2[k]) / fmax(su2[k], 1e-30) < tol))
                converged = 0;
            su2[k] = new_su2;
        }

        double new_s2 = s2 * a_e / fmax(tr_Vinv, 1e-30);
        if (new_s2 < floor_v) new_s2 = floor_v;
        if (!(fabs(new_s2 - s2) / fmax(s2, 1e-30) < tol)) converged = 0;
        s2 = new_s2;

        if (converged) break;
    }

    if (P_out) {
        /* Rebuild V at the fitted components and form the REML projector. */
        for (size_t q = 0; q < (size_t)n * n; q++) P_out[q] = su2[0] * K[0][q];
        for (int k = 1; k < nK; k++)
            for (size_t q = 0; q < (size_t)n * n; q++) P_out[q] += su2[k] * K[k][q];
        for (int i = 0; i < n; i++) P_out[(size_t)i * n + i] += s2;
        int info = 0;
        dpotrf_("L", &n, P_out, &n, &info);
        if (info == 0) dpotri_("L", &n, P_out, &n, &info);
        if (info != 0) { free(yc); free(V); free(Vinv); free(Vy); free(KVy); return -1; }
        symmetrize_full(P_out, n);
        double *v1 = Vy;                   /* reuse: v1 = V^{-1} 1 = row sums */
        double s11 = 0.0;
        for (int i = 0; i < n; i++) {
            const double *Pi = P_out + (size_t)i * n;
            double r = 0.0;
            for (int j = 0; j < n; j++) r += Pi[j];
            v1[i] = r; s11 += r;
        }
        if (!(fabs(s11) > 0.0)) { free(yc); free(V); free(Vinv); free(Vy); free(KVy); return -1; }
        for (int i = 0; i < n; i++)
            for (int j = 0; j < n; j++)
                P_out[(size_t)i * n + j] -= v1[i] * v1[j] / s11;
    }

    /* alpha = A^{-1} y_c with A = sum_k (su2_k/su2_0) K_k + (s2/su2_0) I. */
    double base = fmax(su2[0], 1e-30);
    double lam = s2 / base;
    double *A = V;   /* reuse V buffer */
    memcpy(A, K[0], (size_t)n * n * sizeof(double));
    for (int k = 1; k < nK; k++) {
        double r = su2[k] / base;
        for (size_t q = 0; q < (size_t)n * n; q++) A[q] += r * K[k][q];
    }
    for (int i = 0; i < n; i++) A[(size_t)i * n + i] += lam;
    memcpy(alpha, yc, (size_t)n * sizeof(double));
    int nrhs = 1, info = 0;
    dposv_("L", &n, &nrhs, A, &n, alpha, &n, &info);
    if (info != 0) {
        free(yc); free(V); free(Vinv); free(Vy); free(KVy);
        return -1;
    }

    *out_s2 = s2;
    for (int k = 0; k < nK; k++) out_su2[k] = su2[k];
    if (ratio_out) for (int k = 0; k < nK; k++) ratio_out[k] = su2[k] / base;
    *out_mean_y = mean;
    free(yc); free(V); free(Vinv); free(Vy); free(KVy);
    return 0;
}

/* Single-component form used by every pre-existing arm. */
static int fit_blup(const double *y, const double *K, int n,
                    double *out_s2, double *out_su2,
                    double *alpha, double *out_mean_y)
{
    return fit_blup_m(y, &K, 1, n, out_s2, out_su2, alpha, out_mean_y,
                      NULL, NULL);
}

/* out[c] = (K_cross @ alpha)[c] + mean_y for the (n_cr × n_tr) cross block. */
static void predict_blup(const double *K_cross, int n_cr, int n_tr,
                         const double *alpha, double mean_y, double *out)
{
    for (int c = 0; c < n_cr; c++) {
        const double *Kc = K_cross + (size_t)c * n_tr;
        double s = 0.0;
        for (int j = 0; j < n_tr; j++) s += Kc[j] * alpha[j];
        out[c] = s + mean_y;
    }
}

/* Multi-component form: yhat[c] = sum_k ratio[k] (K_cross[k] @ alpha)[c] + mean_y. */
static void predict_blup_m(const double *const *K_cross, const double *ratio,
                           int nK, int n_cr, int n_tr,
                           const double *alpha, double mean_y, double *out)
{
    for (int c = 0; c < n_cr; c++) out[c] = mean_y;
    for (int k = 0; k < nK; k++) {
        const double *Kk = K_cross[k];
        double r = ratio[k];
        for (int c = 0; c < n_cr; c++) {
            const double *Kc = Kk + (size_t)c * n_tr;
            double s = 0.0;
            for (int j = 0; j < n_tr; j++) s += Kc[j] * alpha[j];
            out[c] += r * s;
        }
    }
}

/* ------------------------------------------------------------------------- *
 * ws² train-tree → full-tree mapping (predict.py:map_train_to_full).
 * For each full tree, find the train tree whose interval contains its midpoint.
 * No pruning by default, so ws2_train[s] is per train tree s directly.
 * ------------------------------------------------------------------------- */
static int *map_ws2_train_to_full(const double *train_breaks, int S_train,
                                  const double *full_breaks, int S_full)
{
    int *idx = malloc((size_t)S_full * sizeof(int));
    if (!idx) return NULL;
    for (int t = 0; t < S_full; t++) {
        double mid = 0.5 * (full_breaks[t] + full_breaks[t + 1]);
        /* searchsorted(train_breaks, mid, side="right") - 1, clamped */
        int lo = 0, hi = S_train + 1;   /* breaks has S_train+1 entries */
        while (lo < hi) {
            int m = (lo + hi) / 2;
            if (train_breaks[m] <= mid) lo = m + 1; else hi = m;
        }
        int k = lo - 1;
        if (k < 0) k = 0;
        if (k > S_train - 1) k = S_train - 1;
        idx[t] = k;
    }
    return idx;
}

/* ------------------------------------------------------------------------- *
 * Kinship block builder.
 *
 * One set of n_train tree passes over the FULL ARG.  Column j = unit vector at
 * train individual j; each pass yields G_s·e_j over all n_total individuals.
 * The visitor accumulates the span-weighted eGRM.  Train rows fill the *_train
 * blocks; test rows the *_cross blocks.  Blocks are row-major.  (It used to
 * accumulate a second, ŵ²-weighted kernel alongside; that fed SARGENT-R alone,
 * retired 2026-08-31.)
 * ------------------------------------------------------------------------- */
typedef struct {
    double *uni;        /* length n_total accumulator for K_uni column */
    int n_total;
} kin_ctx_t;

static void kin_visit(const tree_walk_state_t *st, void *vctx)
{
    kin_ctx_t *c = vctx;
    double span = st->span;
    const double *gv = st->gv;
    /* Only the uniform accumulator remains.  The ŵ²-weighted one (K_sg) fed
     * SARGENT-R alone, retired 2026-08-31; sargent_w / sargent_wf use branch
     * selection instead.  The uniform kernel is independent of ŵ², which is why
     * every surviving arm needs no Stage-2 fit. */
    for (int i = 0; i < c->n_total; i++) c->uni[i] += span * gv[i];
}

static int build_kinship_blocks(const tsk_treeseq_t *full, int n_tr, int n_te,
                                const int *h2i_full, int S_full,
                                double *K_uni_train, double *K_uni_cross,
                                int silent)
{
    int n_total = n_tr + n_te;

    if (!silent)
        fprintf(stderr, "predict: building kinship blocks "
                "(%d tree passes over %d full-ARG trees) ...\n", n_tr, S_full);

    /* Parallelise over COLUMNS, not over samples inside the walker.
     *
     * Each column j is fully independent: it needs only its own probe vector,
     * accumulators and incremental walker state, and it writes a disjoint
     * column of K.  Threading here rather than inside the walker matters
     * because the readout is now a prefix sweep with a serial dependence (a
     * node needs its parent), so it cannot workshare; given a whole column
     * each, threads scale on top of the ~7x work reduction instead of
     * competing with it.  tree_walker_run_inc detects the enclosing team via
     * SARGENT_OMP_IN_PARALLEL() and takes its serial path, so no nesting.
     *
     * Per-thread footprint is small: the inc state is 5 arrays x n_nodes
     * (~240 KB on the production ARG) plus four n_total vectors.
     *
     * OOM or a walker error sets a per-thread flag; `omp for` cannot be
     * broken out of, so remaining iterations fall through cheaply and the
     * failure is reduced across the team. */
    int failed = 0;
    #pragma omp parallel reduction(|:failed)
    {
        tree_cov_matvec_inc_state_t inc;
        memset(&inc, 0, sizeof inc);
        double *in = NULL, *gv = NULL, *uni = NULL;
        int bad = 0;

        if (tree_cov_matvec_inc_state_init(&inc, full, n_total, h2i_full) != 0) {
            bad = 1;
        } else {
            in  = calloc((size_t)n_total, sizeof(double));
            gv  = calloc((size_t)n_total, sizeof(double));
            uni = calloc((size_t)n_total, sizeof(double));
            if (!in || !gv || !uni) bad = 1;
        }

        kin_ctx_t ctx = { .uni = uni, .n_total = n_total };

        #pragma omp for schedule(static)
        for (int j = 0; j < n_tr; j++) {
            if (bad) continue;

            in[j] = 1.0;
            memset(uni, 0, (size_t)n_total * sizeof(double));

            int ret = tree_walker_run_inc(&inc, S_full, /*prune_map=*/NULL,
                                          in, gv, /*pre_visit=*/NULL,
                                          kin_visit, &ctx, /*use_span=*/1);
            in[j] = 0.0;
            if (ret != 0) { bad = 1; continue; }

            for (int i = 0; i < n_tr; i++) {
                K_uni_train[(size_t)i * n_tr + j] = uni[i];
            }
            for (int c = 0; c < n_te; c++) {
                K_uni_cross[(size_t)c * n_tr + j] = uni[n_tr + c];
            }
        }

        free(in); free(gv); free(uni);
        tree_cov_matvec_inc_state_free(&inc);
        failed |= bad;
    }
    if (failed) return -1;

    return 0;
}

/* ------------------------------------------------------------------------- *
 * Genotype loading: X_train (p×n_tr) and X_test (p×n_te), standardized by
 * empirical TRAIN mean and population std (ddof=0), matching predict.py.
 * Aligned to sargent's train-polymorphic variant set (β-indexed).
 * ------------------------------------------------------------------------- */
static int load_genotypes(geno_reader_t *train_gr, geno_reader_t *test_gr,
                          int n_tr, int n_te, int p,
                          double *X_train, double *X_test)
{
    double *buf_tr = malloc((size_t)n_tr * sizeof(double));
    double *buf_te = malloc((size_t)n_te * sizeof(double));
    if (!buf_tr || !buf_te) { free(buf_tr); free(buf_te); return -1; }

    for (int j = 0; j < p; j++) {
        if (geno_seek(train_gr, j) != 0) goto err;
        if (geno_read_variant(train_gr, buf_tr, GENO_RAW) != 1) goto err;

        double mean = 0.0;
        for (int i = 0; i < n_tr; i++) mean += buf_tr[i];
        mean /= n_tr;
        double var = 0.0;
        for (int i = 0; i < n_tr; i++) { double d = buf_tr[i] - mean; var += d * d; }
        var /= n_tr;
        double sd = sqrt(var);
        if (sd == 0.0) sd = 1.0;

        double *Xtr = X_train + (size_t)j * n_tr;
        for (int i = 0; i < n_tr; i++) Xtr[i] = (buf_tr[i] - mean) / sd;

        int file_idx = geno_bed_poly_file_index(train_gr, j);
        if (file_idx < 0) goto err;
        if (geno_bed_read_raw_at(test_gr, file_idx, buf_te) != 0) goto err;
        double *Xte = X_test + (size_t)j * n_te;
        for (int i = 0; i < n_te; i++) Xte[i] = (buf_te[i] - mean) / sd;
    }

    free(buf_tr); free(buf_te);
    return 0;
err:
    free(buf_tr); free(buf_te);
    return -1;
}

/* K_g_train[a,b] = (1/p) Σ_j X_train[j,a] X_train[j,b]  (n_tr×n_tr)
 * K_g_cross[c,b] = (1/p) Σ_j X_test[j,c]  X_train[j,b]  (n_te×n_tr) */
static void snp_grm(const double *X_train, const double *X_test,
                    int n_tr, int n_te, int p,
                    double *K_g_train, double *K_g_cross)
{
    memset(K_g_train, 0, (size_t)n_tr * n_tr * sizeof(double));
    memset(K_g_cross, 0, (size_t)n_te * n_tr * sizeof(double));
    for (int j = 0; j < p; j++) {
        const double *xtr = X_train + (size_t)j * n_tr;
        for (int a = 0; a < n_tr; a++) {
            double xa = xtr[a];
            double *Ka = K_g_train + (size_t)a * n_tr;
            for (int b = 0; b < n_tr; b++) Ka[b] += xa * xtr[b];
        }
        const double *xte = X_test + (size_t)j * n_te;
        for (int c = 0; c < n_te; c++) {
            double xc = xte[c];
            double *Kc = K_g_cross + (size_t)c * n_tr;
            for (int b = 0; b < n_tr; b++) Kc[b] += xc * xtr[b];
        }
    }
    double invp = 1.0 / p;
    for (size_t k = 0; k < (size_t)n_tr * n_tr; k++) K_g_train[k] *= invp;
    for (size_t k = 0; k < (size_t)n_te * n_tr; k++) K_g_cross[k] *= invp;
}

/* Assemble an OUT-OF-FOLD Stage-1 PRS for the training individuals.
 *
 * Fold f holds individuals i with (i % K) == f.  Stage 1 is refitted on the
 * complement and used to score fold f, so no individual is scored by a model
 * that saw them.  Returns 0 on success; on any failure the caller falls back to
 * the in-sample PRS rather than silently producing a half-crossfitted vector.
 *
 * X_train is variant-major and standardized against the full training sample.
 * The subset reader re-standardizes within each fold and publishes the fold's
 * per-variant transform in stats[], which is applied here when scoring the
 * held-out rows -- so fit and prediction share one scale.
 */
static int stage1_oof_prs(sargent_t *sg, const double *X_train,
                          const double *y_train, int p, int n_tr, int K,
                          int silent, double *prs_oof)
{
    int rc = -1;
    int *rows = malloc((size_t)n_tr * sizeof(int));
    double *yfold = malloc((size_t)n_tr * sizeof(double));
    if (!rows || !yfold) { free(rows); free(yfold); return -1; }

    for (int i = 0; i < n_tr; i++) prs_oof[i] = 0.0;

    for (int f = 0; f < K; f++) {
        int nf = 0;
        for (int i = 0; i < n_tr; i++)
            if (i % K != f) { rows[nf] = i; yfold[nf] = y_train[i]; nf++; }
        if (nf < 2) goto done;

        geno_reader_t *gr = geno_open_mem_subset(X_train, p, n_tr, rows, nf);
        if (!gr) goto done;
        if (geno_precompute_stats(gr) != 0) { geno_close(gr); goto done; }

        /* ts is borrowed and only used for the segment map, which Stage 1 does
         * not consult; the fold model inherits the global Gibbs setting. */
        sargent_t *fold = sargent_new(gr, sg->ts, yfold);
        if (!fold) { geno_close(gr); goto done; }

        var_fixed(fold, 1 /* silent: K fits would otherwise flood the log */);

        if (!fold->beta) { sargent_free(fold); geno_close(gr); goto done; }

        /* Score the held-out rows with THIS fold's model, on this fold's scale. */
        for (int i = 0; i < n_tr; i++) {
            if (i % K == f) {
                double acc = 0.0;
                for (int j = 0; j < p; j++) {
                    double sc = gr->stats[j].scale;
                    if (sc == 0.0) continue;          /* monomorphic in-fold */
                    double x = (X_train[(size_t)j * n_tr + i] - gr->stats[j].mean) * sc;
                    acc += fold->beta[j] * x;
                }
                prs_oof[i] = acc;
            }
        }
        sargent_free(fold);
        geno_close(gr);
        if (!silent)
            fprintf(stderr, "stage1-cv: fold %d/%d fitted on %d individuals\n",
                    f + 1, K, nf);
    }
    rc = 0;
done:
    free(rows); free(yfold);
    return rc;
}

/* Out-of-fold ARG-BLUP predictions for the training individuals.
 *
 * The stacking counterpart of stage1_oof_prs: fold f is predicted by a BLUP
 * fitted on the complement, so the weights are learned from predictions no
 * component made on data it saw.  K_uni is subsetted directly -- the fold's own
 * kinship block is K[rows,rows] and the cross block is K[fold,rows].
 */
static int arg_oof_blup(const double *K_uni, const double *y_train,
                        int n_tr, int K, double alpha_pow, double *arg_oof)
{
    int rc = -1;
    int *rows = malloc((size_t)n_tr * sizeof(int));
    int *held = malloc((size_t)n_tr * sizeof(int));
    double *ysub = malloc((size_t)n_tr * sizeof(double));
    double *Ksub = malloc((size_t)n_tr * n_tr * sizeof(double));
    double *Kx   = malloc((size_t)n_tr * n_tr * sizeof(double));
    double *al   = malloc((size_t)n_tr * sizeof(double));
    double *pred = malloc((size_t)n_tr * sizeof(double));
    if (!rows || !held || !ysub || !Ksub || !Kx || !al || !pred) goto done;

    for (int i = 0; i < n_tr; i++) arg_oof[i] = 0.0;

    for (int f = 0; f < K; f++) {
        int nf = 0, nh = 0;
        for (int i = 0; i < n_tr; i++) {
            if (i % K != f) { ysub[nf] = y_train[i]; rows[nf++] = i; }
            else             held[nh++] = i;
        }
        if (nf < 2 || nh < 1) goto done;

        for (int a = 0; a < nf; a++)
            for (int b = 0; b < nf; b++)
                Ksub[(size_t)a * nf + b] = K_uni[(size_t)rows[a] * n_tr + rows[b]];
        for (int c = 0; c < nh; c++)
            for (int b = 0; b < nf; b++)
                Kx[(size_t)c * nf + b] = K_uni[(size_t)held[c] * n_tr + rows[b]];

        double s2f, su2f, meanf;
        if (fit_blup(ysub, Ksub, nf, &s2f, &su2f, al, &meanf) != 0) goto done;
        predict_blup(Kx, nh, nf, al, meanf, pred);
        for (int c = 0; c < nh; c++) arg_oof[held[c]] = pred[c];
    }
    rc = 0;
done:
    free(rows); free(held); free(ysub); free(Ksub); free(Kx); free(al); free(pred);
    (void)alpha_pow;
    return rc;
}

/* Least squares for y ~ 1 + a + b (3x3 normal equations, Gaussian elimination).
 * Returns 0 and fills w[3]; -1 if the system is singular. */
static int stack_weights(const double *y, const double *a, const double *b,
                         int n, double *w)
{
    double M[3][4] = {{0}};
    for (int i = 0; i < n; i++) {
        double v[3] = { 1.0, a[i], b[i] };
        for (int r = 0; r < 3; r++) {
            for (int c = 0; c < 3; c++) M[r][c] += v[r] * v[c];
            M[r][3] += v[r] * y[i];
        }
    }
    for (int c = 0; c < 3; c++) {
        int piv = c;
        for (int r = c + 1; r < 3; r++) if (fabs(M[r][c]) > fabs(M[piv][c])) piv = r;
        if (fabs(M[piv][c]) < 1e-12) return -1;
        if (piv != c) for (int k = 0; k < 4; k++) { double t = M[c][k]; M[c][k] = M[piv][k]; M[piv][k] = t; }
        for (int r = 0; r < 3; r++) {
            if (r == c) continue;
            double fct = M[r][c] / M[c][c];
            for (int k = c; k < 4; k++) M[r][k] -= fct * M[c][k];
        }
    }
    for (int c = 0; c < 3; c++) w[c] = M[c][3] / M[c][c];
    return 0;
}

/* ---- Stage-1 cross-fitting -------------------------------------------------
 *
 * Stage 2 fits the ARG random effect to rtrain = y_train - prs_train.  When
 * prs_train comes from a Stage-1 model fitted on those same individuals, the
 * residual is IN-SAMPLE: a well-fitting Stage 1 absorbs signal and noise alike,
 * the residual looks like white noise, and REML sets the ARG variance to zero.
 * Measured on the phi grid, sargent minus sargent_s was +0.000000 at every
 * level, while sargent_a minus sargent was +0.030 -- the ARG carried real
 * signal that the two-stage arm could not use.  Improving Stage 1 made this
 * WORSE, since a tighter training fit leaves less residual.
 *
 * With --stage1-cv K, prs_train is instead assembled out of fold: Stage 1 is
 * refitted K times, each on K-1 folds, and each individual is scored by a model
 * that never saw them.  prs_TEST is deliberately left as the full-data fit --
 * test individuals were never in training, so there is nothing to cross-fit,
 * and using folds there would only weaken the predictor.  Consequently the
 * arms that never touch rtrain (sargent_s, snp_lmm, sargent_a, sargent_w) are
 * bit-identical with and without the flag.
 * --------------------------------------------------------------------------- */
static int s1_cv_folds = 0;      /* 0 or 1 = off */
static int s1_stack    = 0;      /* --stack: emit the stacked arm */
static int s1_joint    = 0;      /* --joint: beta came from a JOINT fit */

void predict_set_stack(int on) { s1_stack = on ? 1 : 0; }
void predict_set_joint(int on) { s1_joint = on ? 1 : 0; }

void predict_set_stage1_cv(int k) { s1_cv_folds = (k > 1) ? k : 0; }

/* prs[i] = Σ_j X[j,i] beta[j]  (X is p×ncol variant-major; beta length p) */
static void compute_prs(const double *X, const double *beta, int p, int ncol,
                        double *prs)
{
    for (int i = 0; i < ncol; i++) prs[i] = 0.0;
    if (!beta) return;
    for (int j = 0; j < p; j++) {
        const double *Xj = X + (size_t)j * ncol;
        double bj = beta[j];
        if (bj == 0.0) continue;
        for (int i = 0; i < ncol; i++) prs[i] += Xj[i] * bj;
    }
}

/* Simple phenotype reader: n values, last whitespace field per non-# line. */
static double *read_phen_n(const char *path, int n)
{
    FILE *f = fopen(path, "r");
    if (!f) { fprintf(stderr, "predict: cannot open %s\n", path); return NULL; }
    double *y = malloc((size_t)n * sizeof(double));
    if (!y) { fclose(f); return NULL; }
    char line[4096];
    int i = 0;
    while (i < n && fgets(line, sizeof(line), f)) {
        if (line[0] == '#') continue;
        char *tok = strtok(line, " \t\n\r"), *last = NULL;
        while (tok) { last = tok; tok = strtok(NULL, " \t\n\r"); }
        if (!last) continue;
        y[i++] = strtod(last, NULL);
    }
    fclose(f);
    if (i < n) {
        fprintf(stderr, "predict: expected %d phenotypes in %s, got %d\n",
                n, path, i);
        free(y);
        return NULL;
    }
    return y;
}

/* Read the "random" column (index 2) of a true_components tsv (header + rows:
 * id fixed random total). Returns malloc'd array [n] or NULL if absent/error. */
static double *read_true_random(const char *path, int n)
{
    FILE *f = fopen(path, "r");
    if (!f) return NULL;
    double *r = malloc((size_t)n * sizeof(double));
    if (!r) { fclose(f); return NULL; }
    char line[4096];
    int i = 0, first = 1;
    while (fgets(line, sizeof(line), f)) {
        if (first) { first = 0; continue; }        /* header */
        if (line[0] == '#') continue;
        char *tok = strtok(line, " \t\n\r");
        int col = 0; double val = NAN;
        while (tok) { if (col == 2) { val = strtod(tok, NULL); break; }
                      tok = strtok(NULL, " \t\n\r"); col++; }
        if (i < n) r[i++] = val;
    }
    fclose(f);
    if (i < n) { free(r); return NULL; }
    return r;
}

/* ------------------------------------------------------------------------- *
 * Driver.
 * ------------------------------------------------------------------------- */
int sargent_predict(sargent_t *sg,
                    const double *y_train,
                    const char *case_prefix,
                    const char *full_trees_path,
                    const char *test_geno_prefix,
                    const char *test_pheno_path,
                    const char *out_json_path,
                    const branchsel_opts_t *bsel,
                    int silent)
{
    int rc = -1;
    int n_tr = sg->n, p = sg->p;

    /* Stage-2 ws2 is needed ONLY by the full `sargent` arm.  Without it we still
     * report sargent_s, snp_lmm, sargent_a and sargent -- the uniform kinship is
     * independent of ws2 -- and omit `sargent` from the JSON so no consumer can
     * mistake a fallback for the real thing.  Stage 2 is ~88% of a cell, so this
     * is the cheap path (--snp-only --predict) for sweeps that do not need it. */
    const int have_ws2 = (sg->ws2 != NULL);
    if (!have_ws2 && !silent)
        fprintf(stderr, "predict: no Stage-2 ws2; reporting four arms "
                "(sargent_s, snp_lmm, sargent_a, sargent), omitting `sargent`.\n");

    /* ---- load full ARG ---- */
    if (!silent) fprintf(stderr, "predict: loading %s ...\n", full_trees_path);
    tsk_treeseq_t full;
    int ret = tsk_treeseq_load(&full, full_trees_path, 0);
    if (ret < 0) {
        fprintf(stderr, "predict: cannot load %s: %s\n",
                full_trees_path, tsk_strerror(ret));
        return -1;
    }
    int n_total = (int)tsk_treeseq_get_num_individuals(&full);
    int n_te = n_total - n_tr;
    int S_full = (int)tsk_treeseq_get_num_trees(&full);
    int S_train = sg->S_orig;

    if (n_te <= 0) {
        fprintf(stderr, "predict: full ARG has %d individuals but n_train=%d; "
                "need a test set (n_test>0).\n", n_total, n_tr);
        tsk_treeseq_free(&full);
        return -1;
    }
    if (!silent)
        fprintf(stderr, "predict: n_train=%d n_test=%d  trees: train=%d full=%d\n",
                n_tr, n_te, S_train, S_full);

    int *h2i_full = tsk_build_hap_to_indiv(&full);
    if (!h2i_full) { tsk_treeseq_free(&full); return -1; }

    /* ---- map ws2 train-tree → full-tree (only needed for the `sargent` arm) ---- */
    int *ws2_map = NULL;
    double *ws2_full = NULL;
    if (have_ws2) {
        const double *train_breaks = tsk_treeseq_get_breakpoints(sg->ts);
        const double *full_breaks  = tsk_treeseq_get_breakpoints(&full);
        ws2_map = map_ws2_train_to_full(train_breaks, S_train,
                                        full_breaks, S_full);
        ws2_full = malloc((size_t)S_full * sizeof(double));
        if (!ws2_map || !ws2_full) goto cleanup_ts;
        for (int t = 0; t < S_full; t++) ws2_full[t] = sg->ws2[ws2_map[t]];
    }

    /* ---- kinship blocks ---- */
    double *K_uni_train = malloc((size_t)n_tr * n_tr * sizeof(double));
    double *K_uni_cross = malloc((size_t)n_te * n_tr * sizeof(double));
    if (!K_uni_train || !K_uni_cross)
        goto cleanup_kin;

    /* THE expensive step: n_tr tree passes over the full ARG (500 x 1741 here).
     * Skip it outright when no requested arm needs an ARG kinship.  Gating the
     * arms alone was not enough -- with --model snp-lmm this walk still ran and
     * took 7226 s, because skipping the .gall cache removes the shortcut but not
     * the work.  The mask, not the cache, has to decide. */
    const unsigned ARG_ARMS = SG_ARM_SARGENT_A | SG_ARM_SARGENT |
                              SG_ARM_SARGENT_W  | SG_ARM_SARGENT_WF |
                              SG_ARM_STACK;
    int need_arg_kinship = ARM_ON(ARG_ARMS);
    if (need_arg_kinship) {
        if (build_kinship_blocks(&full, n_tr, n_te, h2i_full, S_full,
                                 K_uni_train, K_uni_cross, silent) != 0) {
            fprintf(stderr, "predict: kinship block build failed\n");
            goto cleanup_kin;
        }
    } else if (!silent) {
        fprintf(stderr, "predict: no ARG arm requested; skipping the kinship "
                        "tree walk (%d passes over %d trees)\n", n_tr, S_full);
    }

    /* ---- genotypes ---- */
    if (!silent) fprintf(stderr, "predict: loading genotypes ...\n");
    char test_bed_prefix[4096];
    snprintf(test_bed_prefix, sizeof(test_bed_prefix), "%s", test_geno_prefix);
    geno_reader_t *test_gr = geno_open_bed(test_bed_prefix);
    if (!test_gr) {
        fprintf(stderr, "predict: cannot open test BED %s.bed\n", test_bed_prefix);
        goto cleanup_kin;
    }
    if (test_gr->n_samples != n_te) {
        fprintf(stderr, "predict: test BED has %d samples, expected n_test=%d\n",
                test_gr->n_samples, n_te);
        geno_close(test_gr);
        goto cleanup_kin;
    }

    double *X_train = malloc((size_t)p * n_tr * sizeof(double));
    double *X_test  = malloc((size_t)p * n_te * sizeof(double));
    if (!X_train || !X_test) { geno_close(test_gr); free(X_train); free(X_test); goto cleanup_kin; }
    if (load_genotypes(sg->gr, test_gr, n_tr, n_te, p, X_train, X_test) != 0) {
        fprintf(stderr, "predict: genotype load/alignment failed\n");
        geno_close(test_gr); free(X_train); free(X_test); goto cleanup_kin;
    }
    geno_close(test_gr);

    /* ---- test phenotype and (optional) true random component ---- */
    double *y_test = read_phen_n(test_pheno_path, n_te);
    if (!y_test) { free(X_train); free(X_test); goto cleanup_kin; }
    char tc_path[4096];
    snprintf(tc_path, sizeof(tc_path), "%s.true_components.tsv", test_geno_prefix);
    double *true_random = read_true_random(tc_path, n_te);   /* nullable */

    /* ---- PRS (fixed effect) ---- */
    double *prs_train = malloc((size_t)n_tr * sizeof(double));
    double *prs_test  = malloc((size_t)n_te * sizeof(double));
    /* prs_TEST always comes from the full-data Stage 1: the test individuals
     * were never in training, so there is nothing to cross-fit, and folding it
     * would only weaken the predictor every arm adds back. */
    compute_prs(X_test,  sg->beta, p, n_te, prs_test);

    /* SARGENT-F is a DIFFERENT MODEL (fixed effects only, no ARG term), so it
     * gets its own PRS from its own fit.  Without this it would inherit the
     * joint beta, which is fitted to carry only its share alongside the ARG and
     * so predicts badly alone -- 0.67/0.30/0.20 across phi against 0.78/0.73/0.70
     * for a standalone fit.  Same vector when there is no separate fit. */
    double *prs_test_f = NULL;
    if (sg->beta_f) {
        prs_test_f = malloc((size_t)n_te * sizeof(double));
        if (!prs_test_f) goto cleanup_pred;
        compute_prs(X_test, sg->beta_f, p, n_te, prs_test_f);
    }

    if (s1_cv_folds > 1) {
        if (stage1_oof_prs(sg, X_train, y_train, p, n_tr, s1_cv_folds,
                           silent, prs_train) != 0) {
            fprintf(stderr, "predict: stage1-cv failed; falling back to the "
                            "in-sample PRS\n");
            compute_prs(X_train, sg->beta, p, n_tr, prs_train);
        } else if (!silent) {
            fprintf(stderr, "stage1-cv: %d-fold out-of-fold PRS in use for the "
                            "Stage-2 residual\n", s1_cv_folds);
        }
    } else {
        compute_prs(X_train, sg->beta, p, n_tr, prs_train);
    }
    double mean_y_train = 0.0;
    for (int i = 0; i < n_tr; i++) mean_y_train += y_train[i];
    mean_y_train /= n_tr;

    /* ---- typed-SNP GRM ---- */
    double *K_g_train = malloc((size_t)n_tr * n_tr * sizeof(double));
    double *K_g_cross = malloc((size_t)n_te * n_tr * sizeof(double));
    snp_grm(X_train, X_test, n_tr, n_te, p, K_g_train, K_g_cross);

    /* Working buffers for predictions. */
    double *yhat = malloc((size_t)n_te * sizeof(double));
    double *alpha = malloc((size_t)n_tr * sizeof(double));
    double *rtrain = malloc((size_t)n_tr * sizeof(double));   /* residual y - Xβ */
    for (int i = 0; i < n_tr; i++) rtrain[i] = y_train[i] - prs_train[i];

    metrics_t m_bayes, m_snp, m_arg, m_su, m_w, m_wf, m_stack;
    /* Full-data ARG test prediction, kept for the stacked arm. */
    double *arg_test = NULL;
    double stack_w[3] = {0.0, 0.0, 0.0};
    int    stack_ok = 0;
    double snp_h2 = NAN, arg_h2 = NAN;
    double arg_su2 = NAN, arg_sigma2 = NAN, arg_h2_pve = NAN;
    double su_su2 = NAN, su_sigma2 = NAN;
    double s2, su2, mean_fit;

    /* Branch selection (--branch-select).  P is the REML projector of the
     * uniform-ARG fit in arm (c), which is exactly the background the branch
     * score is taken against; it is filled there and consumed below. */
    const int do_bsel = (bsel != NULL);
    double *P_uni = NULL;
    branchsel_t bs; memset(&bs, 0, sizeof(bs));
    double *K_sel_train = NULL, *K_sel_cross = NULL;
    int    bsel_ok = 0;
    double w_h2 = NAN, w_h2_uni = NAN, w_share_sel = NAN, wf_h2 = NAN;
    if (s1_stack) {
        arg_test = malloc((size_t)n_te * sizeof(double));
        if (!arg_test) goto cleanup_pred;
    }
    if (do_bsel) {
        P_uni = malloc((size_t)n_tr * n_tr * sizeof(double));
        if (!P_uni) goto cleanup_pred;
    }

    /* (a) sargent_s = SARGENT-F: yhat = Xtestᵀβ_F + mean_y, from the
     *     fixed-effects-only fit whenever one was made. */
    {
        const double *pf = sg->beta_f ? prs_test_f : prs_test;
        for (int c = 0; c < n_te; c++) yhat[c] = pf[c] + mean_y_train;
    }
    m_bayes = compute_metrics(yhat, y_test, n_te);

    /* (b) snp_lmm: center K_g, fit on y_train, BLUP.
     *
     * ARM_ON gates the FIT, not merely the JSON write.  That distinction is the
     * whole point of --model: guarding only the emission still pays for the
     * work.  Getting it wrong is not hypothetical -- --model snp-lmm once burned
     * 3.4 hours because the .gall build was skipped while the ARG arm still ran,
     * and without the cached eigenbasis it took the slow uncached path. */
    if (ARM_ON(SG_ARM_SNP_LMM)) {
        center_kinship(K_g_train, n_tr, K_g_cross, n_te);
        if (fit_blup(y_train, K_g_train, n_tr, &s2, &su2, alpha, &mean_fit) != 0) goto cleanup_pred;
        predict_blup(K_g_cross, n_te, n_tr, alpha, mean_fit, yhat);
        m_snp = compute_metrics(yhat, y_test, n_te);
        snp_h2 = su2 / (su2 + s2);
    }

    /* (c) sargent_a: center K_uni, fit on y_train, BLUP.  Keep the RAW K_uni for
     *     the sargent arm, so center a copy. */
    /* SARGENT-W and -WF must be here too: this block's fit_blup_m call is what
     * FILLS P_uni, the REML projector the branch scores are taken against.
     * Gating it out left every score at 0, so selection picked nothing and
     * sargent_w silently fell back to the uniform kinship -- reporting exactly
     * sargent_a's value, which is also the documented signature of k=0 and so
     * looks like a plausible result rather than a failure. */
    if (ARM_ON(SG_ARM_SARGENT_A | SG_ARM_SARGENT | SG_ARM_STACK |
               SG_ARM_SARGENT_W  | SG_ARM_SARGENT_WF)) {
        double *Kt = malloc((size_t)n_tr * n_tr * sizeof(double));
        double *Kc = malloc((size_t)n_te * n_tr * sizeof(double));
        memcpy(Kt, K_uni_train, (size_t)n_tr * n_tr * sizeof(double));
        memcpy(Kc, K_uni_cross, (size_t)n_te * n_tr * sizeof(double));
        center_kinship(Kt, n_tr, Kc, n_te);
        const double *Karr[1] = { Kt };
        if (fit_blup_m(y_train, Karr, 1, n_tr, &s2, &su2, alpha, &mean_fit,
                       NULL, P_uni) != 0) {
            free(Kt); free(Kc); goto cleanup_pred;
        }
        predict_blup(Kc, n_te, n_tr, alpha, mean_fit, yhat);
        m_arg = compute_metrics(yhat, y_test, n_te);
        /* Keep the full-data ARG test prediction: stacking blends it with the
         * PRS rather than fitting one to the other's residual. */
        if (arg_test) memcpy(arg_test, yhat, (size_t)n_te * sizeof(double));
        arg_h2 = su2 / (su2 + s2);
        arg_su2 = su2; arg_sigma2 = s2;
        /* Same realized-variance estimator the joint chain reports, so the two
         * SARGENT arms in the h^2 panel are on one footing.  Kt is the CENTRED
         * train kernel this fit used.  On failure arg_h2 keeps the plug-in. */
        {
            double h2r;
            if (realized_h2_1k(Kt, y_train, n_tr, su2, s2, &h2r) == 0)
                arg_h2_pve = h2r;
        }
        free(Kt); free(Kc);
    }

    /* (d) sargent: residual r; center K_uni (in place now); BLUP; + Xβ.
     * Gated: it consumes K_uni, so leaving it ungated is what forced the kinship
     * walk to run even for --model snp-lmm. */
    if (ARM_ON(SG_ARM_SARGENT | SG_ARM_SARGENT_WF | SG_ARM_STACK)) {
        center_kinship(K_uni_train, n_tr, K_uni_cross, n_te);
        /* JOINT fit: nothing special is needed here, and that is the point.
         *
         * An earlier version passed the joint su2 straight through to the BLUP.
         * That is WRONG on scale: var_joint_gibbs works in the eigenbasis of the
         * .gall G_all, which is neither centred nor normalised, whereas
         * K_uni_train here is double-centred and mean-diag normalised.  The two
         * su2 live in different units (the C printed su2 = 0.00000 against the
         * prototype's 0.955; sigma2 agreed to four decimals, which is what
         * identified the mismatch as scale rather than failure).
         *
         * Re-fitting REML on the residual is both correctly scaled AND sound here,
         * because the residual is now honest: beta came from a joint fit and so did
         * not over-claim the shared signal.  The collapse this whole change exists
         * to fix came from a SEQUENTIAL beta absorbing everything before Stage 2
         * ever saw it -- not from re-fitting per se. */
        if (fit_blup(rtrain, K_uni_train, n_tr, &s2, &su2, alpha, &mean_fit) != 0)
            goto cleanup_pred;
        /* Record the ARG variance components behind sargent.  They were being
         * computed and discarded: --silent suppresses the stderr line, so nothing
         * recorded what was actually fitted.  That is a reproducibility gap on its
         * own -- the Methods describe fitting these jointly -- and they are also
         * what a downstream branch score needs, since P must be built from the same
         * V the score is taken against. */
        su_su2 = su2; su_sigma2 = s2;
        predict_blup(K_uni_cross, n_te, n_tr, alpha, mean_fit, yhat);
        for (int c = 0; c < n_te; c++) yhat[c] += prs_test[c];
        m_su = compute_metrics(yhat, y_test, n_te);
    }

    /* (d2) sargent_stack: BLEND the two predictors instead of fitting one to
     *      the other's residual.
     *
     *      Sequential residualisation asks "what does the ARG add beyond the
     *      PRS?", and the answer is ~0 whenever the two capture the same
     *      heritable signal -- even when the ARG ALONE predicts better.  On the
     *      phi grid sargent equalled sargent_s to six decimals at every level
     *      while sargent_a beat sargent_s by up to 0.05.  Stacking asks the question
     *      that actually matters -- what combination predicts best -- by
     *      regressing y on the two OUT-OF-FOLD training predictions and
     *      applying those weights to the full-data test predictions.  Being a
     *      superset of either component, it cannot do worse than the better one
     *      except through weight-estimation noise. */
    if (s1_stack) {
        double *prs_oof = malloc((size_t)n_tr * sizeof(double));
        double *arg_oof = malloc((size_t)n_tr * sizeof(double));
        if (!prs_oof || !arg_oof) { free(prs_oof); free(arg_oof); goto cleanup_pred; }
        int k = (s1_cv_folds > 1) ? s1_cv_folds : 5;
        if (stage1_oof_prs(sg, X_train, y_train, p, n_tr, k, 1, prs_oof) == 0 &&
            arg_oof_blup(K_uni_train, y_train, n_tr, k, 0.0, arg_oof) == 0 &&
            stack_weights(y_train, prs_oof, arg_oof, n_tr, stack_w) == 0) {
            for (int c = 0; c < n_te; c++)
                yhat[c] = stack_w[0] + stack_w[1] * prs_test[c] + stack_w[2] * arg_test[c];
            m_stack = compute_metrics(yhat, y_test, n_te);
            stack_ok = 1;
            if (!silent)
                fprintf(stderr, "stack: w0=%.4f w_prs=%.4f w_arg=%.4f\n",
                        stack_w[0], stack_w[1], stack_w[2]);
        } else {
            fprintf(stderr, "predict: stacking failed; arm omitted\n");
        }
        free(prs_oof); free(arg_oof);
    }

    /* (e) SARGENT-R retired 2026-08-31.  It was the only consumer of the
     *     Stage-2 learned-ws^2 kinship K_sg, so that accumulator is gone too;
     *     sargent_w / sargent_wf use branch selection instead.  The ws^2 TABLE
     *     is still written (est.w2.tsv) for Blocks 2.1/2.2 -- only this
     *     predictor is retired. */

    /* (f) sargent_w / sargent_wf: branch selection in place of the Stage-2 EM.
     *
     *     Rank branches by the marginal LMM score against the uniform-ARG
     *     background, keep the top k, and fit a two-component REML
     *     (K_uni + K_sel).  sargent_w is the pure-ARG arm — the direct
     *     counterpart of sargent_a, which it must reproduce at k = 0 since the
     *     second component simply disappears.  sargent_wf adds the Stage-1 PRS,
     *     making it the counterpart of sargent / sargent.
     *
     *     K_uni_train / K_uni_cross are already centered in place by (d). */
    if (do_bsel) {
        branchsel_opts_t o = *bsel;
        o.silent = silent;
        /* Validation hook: dump the projector and the centered uniform kinship
         * it was built from, so an external reference can re-derive P itself
         * rather than taking it on trust.  Off unless the env var is set. */
        {
            const char *mp = getenv("SARGENT_BRANCHSEL_DUMP_MAT");
            if (mp && *mp) {
                char pth[4200];
                snprintf(pth, sizeof(pth), "%s.P.bin", mp);
                FILE *pf = fopen(pth, "wb");
                if (pf) { fwrite(P_uni, sizeof(double), (size_t)n_tr * n_tr, pf); fclose(pf); }
                snprintf(pth, sizeof(pth), "%s.Kuni.bin", mp);
                FILE *kf = fopen(pth, "wb");
                if (kf) { fwrite(K_uni_train, sizeof(double), (size_t)n_tr * n_tr, kf); fclose(kf); }
                snprintf(pth, sizeof(pth), "%s.Kuni_cross.bin", mp);
                FILE *cf = fopen(pth, "wb");
                if (cf) { fwrite(K_uni_cross, sizeof(double), (size_t)n_te * n_tr, cf); fclose(cf); }
                snprintf(pth, sizeof(pth), "%s.meta.json", mp);
                FILE *mf = fopen(pth, "w");
                if (mf) {
                    fprintf(mf, "{\"n_train\": %d, \"n_total\": %d, \"su2\": ", n_tr, n_total);
                    json_num(mf, arg_su2);
                    fprintf(mf, ", \"sigma2\": "); json_num(mf, arg_sigma2);
                    fprintf(mf, "}\n");
                    fclose(mf);
                }
                if (!silent)
                    fprintf(stderr, "branchsel: wrote projector/kinship dumps to %s.*\n", mp);
            }
        }
        double *y_c = malloc((size_t)n_tr * sizeof(double));
        if (!y_c) goto cleanup_pred;
        for (int i = 0; i < n_tr; i++) y_c[i] = y_train[i] - mean_y_train;
        int brc = branchsel_run(&full, h2i_full, n_total, n_tr, P_uni, y_c,
                                &o, &bs);
        free(y_c);
        if (brc != 0) {
            fprintf(stderr, "predict: branch selection failed\n");
            goto cleanup_pred;
        }
        /* k = 0 is the no-selection control: the second component simply is not
         * there, so the fit is the SAME fit_blup_m call on the SAME centered
         * K_uni that arm (c) makes.  sargent_w then equals sargent_a by
         * construction, which is what makes the arm no-regret. */
        /* Optional dump of the SELECTED-BRANCH carrier matrix Z.
         *
         * K_sel = Z Z' / k is rank k, so a joint (rather than sequential) fit
         * can treat the selected branches as k extra design columns with a
         * Gaussian prior instead of forming the kernel -- which keeps the
         * eigenbasis of K_uniform diagonal and makes the joint fit cheap.  The
         * Python prototype needs Z to test that; the existing
         * SARGENT_BRANCHSEL_DUMP writes only scores and fingerprints.
         * Columns are standardized by TRAIN mean and population sd, matching
         * what branchsel_kinship consumes. */
        {
            const char *zp = getenv("SARGENT_BRANCHSEL_DUMP_Z");
            if (zp && *zp && bs.k > 0 && bs.Z) {
                char pth[4200];
                snprintf(pth, sizeof(pth), "%s.Z.bin", zp);
                FILE *zf = fopen(pth, "wb");
                if (zf) {
                    fwrite(bs.Z, sizeof(double),
                           (size_t)bs.n_total * (size_t)bs.k, zf);
                    fclose(zf);
                }
                snprintf(pth, sizeof(pth), "%s.Z.meta.json", zp);
                FILE *mf = fopen(pth, "w");
                if (mf) {
                    fprintf(mf, "{\"n_total\": %d, \"n_train\": %d, \"k\": %d,"
                                " \"layout\": \"row-major n_total x k\"}\n",
                            bs.n_total, bs.n_train, bs.k);
                    fclose(mf);
                }
                if (!silent)
                    fprintf(stderr, "branchsel: wrote Z (%d x %d) to %s.Z.bin\n",
                            bs.n_total, bs.k, zp);
            }
        }

        int nK = (bs.k > 0) ? 2 : 1;
        if (nK == 2) {
            K_sel_train = malloc((size_t)n_tr * n_tr * sizeof(double));
            K_sel_cross = malloc((size_t)n_te * n_tr * sizeof(double));
            if (!K_sel_train || !K_sel_cross) goto cleanup_pred;
            branchsel_kinship(&bs, K_sel_train, K_sel_cross);
            center_kinship(K_sel_train, n_tr, K_sel_cross, n_te);
        }

        const double *Ktr2[2] = { K_uni_train, K_sel_train };
        const double *Kcr2[2] = { K_uni_cross, K_sel_cross };
        double su2v[2] = { 0.0, 0.0 }, ratio[2] = { 0.0, 0.0 };

        if (fit_blup_m(y_train, Ktr2, nK, n_tr, &s2, su2v, alpha, &mean_fit,
                       ratio, NULL) != 0) goto cleanup_pred;
        predict_blup_m(Kcr2, ratio, nK, n_te, n_tr, alpha, mean_fit, yhat);
        m_w = compute_metrics(yhat, y_test, n_te);
        /* The selected component's variance is upward-biased by selection, so
         * w_h2_uni (the uniform component alone) is the defensible estimate and
         * w_h2 is reported only for comparability with the prototype. */
        w_h2        = (su2v[0] + su2v[1]) / (su2v[0] + su2v[1] + s2);
        w_h2_uni    = su2v[0] / (su2v[0] + s2);
        w_share_sel = su2v[1] / fmax(su2v[0] + su2v[1], 1e-300);

        if (fit_blup_m(rtrain, Ktr2, nK, n_tr, &s2, su2v, alpha, &mean_fit,
                       ratio, NULL) != 0) goto cleanup_pred;
        predict_blup_m(Kcr2, ratio, nK, n_te, n_tr, alpha, mean_fit, yhat);
        for (int c = 0; c < n_te; c++) yhat[c] += prs_test[c];
        m_wf = compute_metrics(yhat, y_test, n_te);
        wf_h2 = (su2v[0] + su2v[1]) / (su2v[0] + su2v[1] + s2);
        bsel_ok = 1;
    }

    /* ---- write JSON ---- */
    {
        FILE *f = fopen(out_json_path, "w");
        if (!f) { fprintf(stderr, "predict: cannot open %s\n", out_json_path); goto cleanup_pred; }
        fprintf(f, "{\n");
        fprintf(f, "  \"case\": \"%s\",\n", case_prefix);
        fprintf(f, "  \"n_train\": %d,\n", n_tr);
        fprintf(f, "  \"n_test\": %d,\n", n_te);
        fprintf(f, "  \"n_trees_full\": %d,\n", S_full);
        fprintf(f, "  \"n_trees_train\": %d,\n", S_train);
        /* ARM_IF emits only when the arm was actually computed.  Without it,
         * unrequested arms were written as FABRICATED ZEROS -- --model snp-lmm
         * produced sargent_a = 0 and sargent = 0 beside a correct snp_lmm, and
         * an aggregator reading those keys sees a real, catastrophic result
         * rather than missing data.  Same failure as the `sargent` arm. */
        #define ARM_IF(bit, name, mm) do { if (ARM_ON(bit)) { ARM(name, mm); } } while (0)
        #define ARM(name, mm) do { \
            fprintf(f, "  \"" name "\": {\n"); \
            fprintf(f, "    \"r\": ");        json_num(f, (mm).r);        fprintf(f, ",\n"); \
            fprintf(f, "    \"r2\": ");       json_num(f, (mm).r2);       fprintf(f, ",\n"); \
            fprintf(f, "    \"r2_calib\": "); json_num(f, (mm).r2_calib); fprintf(f, ",\n"); \
            fprintf(f, "    \"mse\": ");      json_num(f, (mm).mse);      fprintf(f, "\n"); \
            fprintf(f, "  },\n"); \
        } while (0)
        ARM_IF(SG_ARM_SARGENT_S, "sargent_s", m_bayes);
        ARM_IF(SG_ARM_SNP_LMM, "snp_lmm", m_snp);
        fprintf(f, "  \"snp_lmm_h2\": "); json_num(f, snp_h2); fprintf(f, ",\n");
        ARM_IF(SG_ARM_SARGENT_A, "sargent_a", m_arg);
        if (arg_h2_pve == arg_h2_pve) {
            fprintf(f, "  \"sargent_a_pve\": "); json_num(f, arg_h2_pve); fprintf(f, ",\n");
        }
        fprintf(f, "  \"sargent_a_h2\": "); json_num(f, arg_h2); fprintf(f, ",\n");
        /* The variance components of this fit define V, hence the REML projector
         * P that branch scores are taken against.  Emitted so an external
         * reference implementation can rebuild exactly the same background. */
        fprintf(f, "  \"sargent_a_su2\": "); json_num(f, arg_su2); fprintf(f, ",\n");
        fprintf(f, "  \"sargent_a_sigma2\": "); json_num(f, arg_sigma2); fprintf(f, ",\n");
        ARM_IF(SG_ARM_SARGENT, "sargent", m_su);
        /* BSLMM-style PVE from the joint chain -- a realized-variance estimate,
         * emitted beside the plug-in so the two can be compared directly on the
         * same run.  Decomposed into its fixed and random parts. */
        if (sg->joint_pve > 0.0) {
            fprintf(f, "  \"sargent_pve\": ");        json_num(f, sg->joint_pve);
            fprintf(f, ",\n  \"sargent_pve_fixed\": ");  json_num(f, sg->joint_pve_fix);
            fprintf(f, ",\n  \"sargent_pve_random\": "); json_num(f, sg->joint_pve_rand);
            fprintf(f, ",\n");
        }
        fprintf(f, "  \"sargent_su2\": ");    json_num(f, su_su2);
        fprintf(f, ",\n  \"sargent_sigma2\": "); json_num(f, su_sigma2);
        fprintf(f, ",\n");
        if (stack_ok) {
            ARM("sargent_stack", m_stack);
            fprintf(f, "  \"sargent_stack_w_prs\": "); json_num(f, stack_w[1]);
            fprintf(f, ",\n  \"sargent_stack_w_arg\": "); json_num(f, stack_w[2]);
            fprintf(f, ",\n");
        }
        /* `sargent` is OMITTED (not zeroed, not aliased to sargent) when there
         * is no Stage-2 ws2, so a consumer that expects it fails loudly rather
         * than silently reading a stand-in.  A flag records why. */
        /* have_ws2 alone is not enough: with --model excluding SARGENT-R the ws2
         * buffer can exist while Stage 2 never ran, and the arm was emitted as a
         * FABRICATED ZERO (r2 = 0, r2_calib = 0, sargent_sigma2 = nan).  An
         * aggregator reading that key sees 0.0 and treats it as a real, terrible
         * result -- it does not look like missing data.  Emit the arm only when
         * it was actually computed; otherwise take the existing omission path,
         * which downstream already understands (aggregate_fig23.py checks for
         * sargent_omitted). */
        fprintf(f, "  \"sargent_omitted\": \"SARGENT-R retired 2026-08-31; "
                   "arms reported: sargent_s, snp_lmm, sargent_a, sargent\"\n");
        if (bsel_ok) {
            fprintf(f, ",\n");
            ARM("sargent_w", m_w);
            ARM("sargent_wf", m_wf);
            fprintf(f, "  \"sargent_w_h2\": ");        json_num(f, w_h2);
            fprintf(f, ",\n  \"sargent_w_h2_uniform\": "); json_num(f, w_h2_uni);
            fprintf(f, ",\n  \"sargent_w_share_sel\": ");  json_num(f, w_share_sel);
            fprintf(f, ",\n  \"sargent_wf_h2\": ");    json_num(f, wf_h2);
            fprintf(f, ",\n  \"sargent_w_k\": %d", bs.k);
            fprintf(f, ",\n  \"sargent_w_n_branch\": %d", bs.n_branch);
            fprintf(f, ",\n  \"sargent_w_n_trees_scanned\": %d", bs.n_trees);
            fprintf(f, ",\n  \"sargent_w_score_top\": ");
            json_num(f, bs.k > 0 ? bs.score[0] : NAN);
            fprintf(f, ",\n  \"sargent_w_score_min\": ");
            json_num(f, bs.k > 0 ? bs.score[bs.k - 1] : NAN);
            fprintf(f, ",\n  \"sargent_w_seconds\": ");    json_num(f, bs.seconds);
            fprintf(f, "\n");
        }
        #undef ARM
        fprintf(f, "}\n");
        fclose(f);
        if (!silent) fprintf(stderr, "predict: wrote %s\n", out_json_path);
    }

    if (!silent) {
        fprintf(stderr, "\n=== Out-of-sample prediction (native) ===\n");
        fprintf(stderr, "  sargent_s     R2=%+.3f r=%+.3f\n", m_bayes.r2, m_bayes.r);
        fprintf(stderr, "  snp_lmm    R2=%+.3f r=%+.3f\n", m_snp.r2, m_snp.r);
        fprintf(stderr, "  sargent_a    R2=%+.3f r=%+.3f\n", m_arg.r2, m_arg.r);
        fprintf(stderr, "  sargent  R2=%+.3f r=%+.3f\n", m_su.r2, m_su.r);
        if (bsel_ok) {
            fprintf(stderr, "  sargent_w  R2=%+.3f r=%+.3f   "
                    "(k=%d of %d branches; dr2_calib vs sargent_a = %+.4f)\n",
                    m_w.r2, m_w.r, bs.k, bs.n_branch,
                    m_w.r2_calib - m_arg.r2_calib);
            fprintf(stderr, "  sargent_wf R2=%+.3f r=%+.3f\n", m_wf.r2, m_wf.r);
        }
    }
    rc = 0;

cleanup_pred:
    branchsel_free(&bs);
    free(P_uni); free(K_sel_train); free(K_sel_cross);
    free(yhat); free(alpha); free(rtrain);
    free(prs_train); free(prs_test); free(prs_test_f);
    free(K_g_train); free(K_g_cross);
    free(X_train); free(X_test);
    free(y_test); free(true_random);
cleanup_kin:
    free(K_uni_train); free(K_uni_cross);
cleanup_ts:
    free(ws2_map); free(ws2_full);
    free(h2i_full);
    tsk_treeseq_free(&full);
    return rc;
}
