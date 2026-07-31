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

#include <tskit.h>
#include "predict.h"
#include "geno.h"
#include "sargent.h"

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
 * fit_blup (predict.py:fit_blup): REML-EM variance components + BLUP alpha.
 * K is row-major n×n (centered).  Returns 0 on success.
 * alpha[n] receives (K + lam I)^{-1} (y - mean_y); *mean_y = mean(y).
 * ------------------------------------------------------------------------- */
static int fit_blup(const double *y, const double *K, int n,
                    double *out_s2, double *out_su2,
                    double *alpha, double *out_mean_y)
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

    double su2 = 0.5 * var_y, s2 = 0.5 * var_y;

    double *V     = malloc((size_t)n * n * sizeof(double));
    double *Vinv  = malloc((size_t)n * n * sizeof(double));
    double *Vy    = malloc((size_t)n * sizeof(double));
    double *KVy   = malloc((size_t)n * sizeof(double));
    if (!yc || !V || !Vinv || !Vy || !KVy) {
        free(yc); free(V); free(Vinv); free(Vy); free(KVy);
        return -1;
    }

    for (int iter = 0; iter < max_iter; iter++) {
        /* V = su2*K + s2*I */
        for (size_t k = 0; k < (size_t)n * n; k++) V[k] = su2 * K[k];
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
        matvec(K, Vy, KVy, n);             /* K Vinv_y */

        double a_u = 0.0, a_e = 0.0;
        for (int i = 0; i < n; i++) { a_u += KVy[i] * Vy[i]; a_e += Vy[i] * Vy[i]; }

        double tr_VinvK = 0.0, tr_Vinv = 0.0;
        for (int i = 0; i < n; i++) {
            const double *Vi = Vinv + (size_t)i * n;
            const double *Ki = K + (size_t)i * n;
            double row = 0.0;
            for (int j = 0; j < n; j++) row += Vi[j] * Ki[j];
            tr_VinvK += row;
            tr_Vinv += Vi[i];
        }

        double new_su2 = su2 * a_u / fmax(tr_VinvK, 1e-30);
        double new_s2  = s2  * a_e / fmax(tr_Vinv,  1e-30);
        double floor_v = 1e-8 * var_y;
        if (new_su2 < floor_v) new_su2 = floor_v;
        if (new_s2  < floor_v) new_s2  = floor_v;

        int converged =
            fabs(new_su2 - su2) / fmax(su2, 1e-30) < tol &&
            fabs(new_s2  - s2 ) / fmax(s2,  1e-30) < tol;
        su2 = new_su2; s2 = new_s2;
        if (converged) break;
    }

    /* alpha = (K + lam I)^{-1} y_c, lam = s2/su2 */
    double lam = s2 / fmax(su2, 1e-30);
    double *A = V;   /* reuse V buffer */
    memcpy(A, K, (size_t)n * n * sizeof(double));
    for (int i = 0; i < n; i++) A[(size_t)i * n + i] += lam;
    memcpy(alpha, yc, (size_t)n * sizeof(double));
    int nrhs = 1, info = 0;
    dposv_("L", &n, &nrhs, A, &n, alpha, &n, &info);
    if (info != 0) {
        free(yc); free(V); free(Vinv); free(Vy); free(KVy);
        return -1;
    }

    *out_s2 = s2; *out_su2 = su2; *out_mean_y = mean;
    free(yc); free(V); free(Vinv); free(Vy); free(KVy);
    return 0;
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
 * The visitor accumulates two weightings simultaneously (span-weighted eGRM and
 * learned ŵ²·span).  Train rows fill the *_train blocks; test rows the *_cross
 * blocks.  Blocks are row-major.
 * ------------------------------------------------------------------------- */
typedef struct {
    double *uni;        /* length n_total accumulator for K_uni column */
    double *sg;         /* length n_total accumulator for K_sg column; NULL when
                         * ws2 is absent (--snp-only): the learned arm is skipped */
    const double *ws2_full;   /* NULL iff sg is NULL */
    int n_total;
} kin_ctx_t;

static void kin_visit(const tree_walk_state_t *st, void *vctx)
{
    kin_ctx_t *c = vctx;
    double span = st->span;
    const double *gv = st->gv;
    /* The uniform accumulator is independent of ws2, which is why the four arms
     * that use it (bayesR, snp_lmm, arg_lmm, sargent_u) need no Stage-2 fit. */
    if (!c->sg) {
        for (int i = 0; i < c->n_total; i++) c->uni[i] += span * gv[i];
        return;
    }
    double w = c->ws2_full[st->s_orig] * span;
    for (int i = 0; i < c->n_total; i++) {
        double g = gv[i];
        c->uni[i] += span * g;
        c->sg[i]  += w * g;
    }
}

static int build_kinship_blocks(const tsk_treeseq_t *full, int n_tr, int n_te,
                                const int *h2i_full, int S_full,
                                const double *ws2_full,
                                double *K_uni_train, double *K_uni_cross,
                                double *K_sg_train,  double *K_sg_cross,
                                int silent)
{
    int n_total = n_tr + n_te;
    tree_cov_matvec_inc_state_t inc;
    if (tree_cov_matvec_inc_state_init(&inc, full, n_total, h2i_full) != 0)
        return -1;

    /* K_sg_train/K_sg_cross are NULL when there is no Stage-2 ws2 fit; then the
     * learned accumulator is not allocated and the walker skips it entirely. */
    int have_sg = (ws2_full != NULL && K_sg_train != NULL && K_sg_cross != NULL);

    double *in  = calloc((size_t)n_total, sizeof(double));
    double *gv  = calloc((size_t)n_total, sizeof(double));
    double *uni = calloc((size_t)n_total, sizeof(double));
    double *sgc = have_sg ? calloc((size_t)n_total, sizeof(double)) : NULL;
    if (!in || !gv || !uni || (have_sg && !sgc)) {
        free(in); free(gv); free(uni); free(sgc);
        tree_cov_matvec_inc_state_free(&inc);
        return -1;
    }

    kin_ctx_t ctx = { .uni = uni, .sg = sgc,
                      .ws2_full = have_sg ? ws2_full : NULL,
                      .n_total = n_total };

    if (!silent)
        fprintf(stderr, "predict: building kinship blocks "
                "(%d tree passes over %d full-ARG trees) ...\n", n_tr, S_full);

    for (int j = 0; j < n_tr; j++) {
        if (j > 0) in[j - 1] = 0.0;
        in[j] = 1.0;
        memset(uni, 0, (size_t)n_total * sizeof(double));
        if (have_sg) memset(sgc, 0, (size_t)n_total * sizeof(double));

        int ret = tree_walker_run_inc(&inc, S_full, /*prune_map=*/NULL,
                                      in, gv, /*pre_visit=*/NULL,
                                      kin_visit, &ctx, /*use_span=*/1);
        if (ret != 0) {
            free(in); free(gv); free(uni); free(sgc);
            tree_cov_matvec_inc_state_free(&inc);
            return -1;
        }

        for (int i = 0; i < n_tr; i++) {
            K_uni_train[(size_t)i * n_tr + j] = uni[i];
            if (have_sg) K_sg_train[(size_t)i * n_tr + j] = sgc[i];
        }
        for (int c = 0; c < n_te; c++) {
            K_uni_cross[(size_t)c * n_tr + j] = uni[n_tr + c];
            if (have_sg) K_sg_cross[(size_t)c * n_tr + j] = sgc[n_tr + c];
        }
    }

    free(in); free(gv); free(uni); free(sgc);
    tree_cov_matvec_inc_state_free(&inc);
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
                    int silent)
{
    int rc = -1;
    int n_tr = sg->n, p = sg->p;

    /* Stage-2 ws2 is needed ONLY by the full `sargent` arm.  Without it we still
     * report bayesR, snp_lmm, arg_lmm and sargent_u -- the uniform kinship is
     * independent of ws2 -- and omit `sargent` from the JSON so no consumer can
     * mistake a fallback for the real thing.  Stage 2 is ~88% of a cell, so this
     * is the cheap path (--snp-only --predict) for sweeps that do not need it. */
    const int have_ws2 = (sg->ws2 != NULL);
    if (!have_ws2 && !silent)
        fprintf(stderr, "predict: no Stage-2 ws2; reporting four arms "
                "(bayesR, snp_lmm, arg_lmm, sargent_u), omitting `sargent`.\n");

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
    double *K_sg_train  = have_ws2 ? malloc((size_t)n_tr * n_tr * sizeof(double)) : NULL;
    double *K_sg_cross  = have_ws2 ? malloc((size_t)n_te * n_tr * sizeof(double)) : NULL;
    if (!K_uni_train || !K_uni_cross || (have_ws2 && (!K_sg_train || !K_sg_cross)))
        goto cleanup_kin;

    if (build_kinship_blocks(&full, n_tr, n_te, h2i_full, S_full, ws2_full,
                             K_uni_train, K_uni_cross,
                             K_sg_train, K_sg_cross, silent) != 0) {
        fprintf(stderr, "predict: kinship block build failed\n");
        goto cleanup_kin;
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
    compute_prs(X_train, sg->beta, p, n_tr, prs_train);
    compute_prs(X_test,  sg->beta, p, n_te, prs_test);
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

    metrics_t m_bayes, m_snp, m_arg, m_su, m_sg;
    double snp_h2 = NAN, arg_h2 = NAN, sargent_sigma2 = NAN, sargent_rand_r = NAN;
    double s2, su2, mean_fit;

    /* (a) bayesR: yhat = Xtestᵀβ + mean_y */
    for (int c = 0; c < n_te; c++) yhat[c] = prs_test[c] + mean_y_train;
    m_bayes = compute_metrics(yhat, y_test, n_te);

    /* (b) snp_lmm: center K_g, fit on y_train, BLUP */
    center_kinship(K_g_train, n_tr, K_g_cross, n_te);
    if (fit_blup(y_train, K_g_train, n_tr, &s2, &su2, alpha, &mean_fit) != 0) goto cleanup_pred;
    predict_blup(K_g_cross, n_te, n_tr, alpha, mean_fit, yhat);
    m_snp = compute_metrics(yhat, y_test, n_te);
    snp_h2 = su2 / (su2 + s2);

    /* (c) arg_lmm: center K_uni, fit on y_train, BLUP.  Keep the RAW K_uni for
     *     the sargent_u arm, so center a copy. */
    {
        double *Kt = malloc((size_t)n_tr * n_tr * sizeof(double));
        double *Kc = malloc((size_t)n_te * n_tr * sizeof(double));
        memcpy(Kt, K_uni_train, (size_t)n_tr * n_tr * sizeof(double));
        memcpy(Kc, K_uni_cross, (size_t)n_te * n_tr * sizeof(double));
        center_kinship(Kt, n_tr, Kc, n_te);
        if (fit_blup(y_train, Kt, n_tr, &s2, &su2, alpha, &mean_fit) != 0) {
            free(Kt); free(Kc); goto cleanup_pred;
        }
        predict_blup(Kc, n_te, n_tr, alpha, mean_fit, yhat);
        m_arg = compute_metrics(yhat, y_test, n_te);
        arg_h2 = su2 / (su2 + s2);
        free(Kt); free(Kc);
    }

    /* (d) sargent_u: residual r; center K_uni (in place now); BLUP; + Xβ */
    center_kinship(K_uni_train, n_tr, K_uni_cross, n_te);
    if (fit_blup(rtrain, K_uni_train, n_tr, &s2, &su2, alpha, &mean_fit) != 0) goto cleanup_pred;
    predict_blup(K_uni_cross, n_te, n_tr, alpha, mean_fit, yhat);
    for (int c = 0; c < n_te; c++) yhat[c] += prs_test[c];
    m_su = compute_metrics(yhat, y_test, n_te);

    /* (e) sargent: residual r; center K_sg; BLUP; + Xβ.
     *     Skipped entirely without a Stage-2 ws2 fit -- this is the only arm
     *     that needs it, and the JSON below omits it rather than substituting. */
    if (have_ws2) {
    center_kinship(K_sg_train, n_tr, K_sg_cross, n_te);
    if (fit_blup(rtrain, K_sg_train, n_tr, &s2, &su2, alpha, &mean_fit) != 0) goto cleanup_pred;
    predict_blup(K_sg_cross, n_te, n_tr, alpha, mean_fit, yhat);
    for (int c = 0; c < n_te; c++) yhat[c] += prs_test[c];
    m_sg = compute_metrics(yhat, y_test, n_te);
    sargent_sigma2 = s2;
    if (true_random) {
        /* corr(yhat - Xtestβ, true_random) = corr(K_sg_cross·α + mean_r, ...) */
        double *uhat = malloc((size_t)n_te * sizeof(double));
        for (int c = 0; c < n_te; c++) uhat[c] = yhat[c] - prs_test[c];
        double mu = 0.0, mr = 0.0;
        for (int c = 0; c < n_te; c++) { mu += uhat[c]; mr += true_random[c]; }
        mu /= n_te; mr /= n_te;
        double sxy = 0.0, sxx = 0.0, syy = 0.0;
        for (int c = 0; c < n_te; c++) {
            double a = uhat[c] - mu, b = true_random[c] - mr;
            sxy += a * b; sxx += a * a; syy += b * b;
        }
        sargent_rand_r = (sxx > 0 && syy > 0) ? sxy / sqrt(sxx * syy) : NAN;
        free(uhat);
    }
    }   /* end if (have_ws2) -- the `sargent` arm */

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
        #define ARM(name, mm) do { \
            fprintf(f, "  \"" name "\": {\n"); \
            fprintf(f, "    \"r\": ");        json_num(f, (mm).r);        fprintf(f, ",\n"); \
            fprintf(f, "    \"r2\": ");       json_num(f, (mm).r2);       fprintf(f, ",\n"); \
            fprintf(f, "    \"r2_calib\": "); json_num(f, (mm).r2_calib); fprintf(f, ",\n"); \
            fprintf(f, "    \"mse\": ");      json_num(f, (mm).mse);      fprintf(f, "\n"); \
            fprintf(f, "  },\n"); \
        } while (0)
        ARM("bayesR", m_bayes);
        ARM("snp_lmm", m_snp);
        fprintf(f, "  \"snp_lmm_h2\": "); json_num(f, snp_h2); fprintf(f, ",\n");
        ARM("arg_lmm", m_arg);
        fprintf(f, "  \"arg_lmm_h2\": "); json_num(f, arg_h2); fprintf(f, ",\n");
        ARM("sargent_u", m_su);
        /* `sargent` is OMITTED (not zeroed, not aliased to sargent_u) when there
         * is no Stage-2 ws2, so a consumer that expects it fails loudly rather
         * than silently reading a stand-in.  A flag records why. */
        if (have_ws2) {
            ARM("sargent", m_sg);
            fprintf(f, "  \"sargent_sigma2\": "); json_num(f, sargent_sigma2);
            if (true_random) {
                fprintf(f, ",\n  \"sargent_rand_r\": "); json_num(f, sargent_rand_r);
            }
            fprintf(f, "\n");
        } else {
            fprintf(f, "  \"sargent_omitted\": \"no Stage-2 ws2 fit "
                       "(--snp-only); arms reported: bayesR, snp_lmm, arg_lmm, "
                       "sargent_u\"\n");
        }
        #undef ARM
        fprintf(f, "}\n");
        fclose(f);
        if (!silent) fprintf(stderr, "predict: wrote %s\n", out_json_path);
    }

    if (!silent) {
        fprintf(stderr, "\n=== Out-of-sample prediction (native) ===\n");
        fprintf(stderr, "  bayesR     R2=%+.3f r=%+.3f\n", m_bayes.r2, m_bayes.r);
        fprintf(stderr, "  snp_lmm    R2=%+.3f r=%+.3f\n", m_snp.r2, m_snp.r);
        fprintf(stderr, "  arg_lmm    R2=%+.3f r=%+.3f\n", m_arg.r2, m_arg.r);
        fprintf(stderr, "  sargent_u  R2=%+.3f r=%+.3f\n", m_su.r2, m_su.r);
        if (have_ws2)
            fprintf(stderr, "  sargent    R2=%+.3f r=%+.3f\n", m_sg.r2, m_sg.r);
        else
            fprintf(stderr, "  sargent    (omitted: no Stage-2 ws2)\n");
    }
    rc = 0;

cleanup_pred:
    free(yhat); free(alpha); free(rtrain);
    free(prs_train); free(prs_test);
    free(K_g_train); free(K_g_cross);
    free(X_train); free(X_test);
    free(y_test); free(true_random);
cleanup_kin:
    free(K_uni_train); free(K_uni_cross);
    free(K_sg_train); free(K_sg_cross);
cleanup_ts:
    free(ws2_map); free(ws2_full);
    free(h2i_full);
    tsk_treeseq_free(&full);
    return rc;
}
