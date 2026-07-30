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
#include <stdint.h>

#include "gall.h"
#include "geno.h"    /* tree_cov_matvec_inc_state_t, tree_walker_run_inc */

/* -------------------------------------------------------------------------
 * Portable LAPACK dsyev declaration.
 * Signature is identical on macOS/Accelerate, Linux/OpenBLAS, LAPACK ref.
 * ---------------------------------------------------------------------- */
extern void dsyev_(const char *jobz, const char *uplo,
                   const int *n, double *a, const int *lda,
                   double *w, double *work, const int *lwork, int *info);

/* -------------------------------------------------------------------------
 * Visitor used to accumulate G_all[:,j] = Σ_s G_s e_j during one tree walk.
 * ---------------------------------------------------------------------- */
typedef struct {
    double *G;    /* n×n column-major accumulator */
    int     n;
    int     col;  /* which column (= which unit-vector basis index) */
} gall_acc_ctx_t;

static void gall_acc_visit(const tree_walk_state_t *st, void *vctx)
{
    gall_acc_ctx_t *c = vctx;
    double *col = c->G + (size_t)c->col * c->n;
    const double *gv = st->gv;
    double span = st->span;            /* span-weighted eGRM: G_all = Sum_s span_s G_s */
    int n = c->n;
    for (int i = 0; i < n; i++)
        col[i] += span * gv[i];
}

/* -------------------------------------------------------------------------
 * gall_build: n incremental tree-sequence passes then LAPACK eigendecompose.
 * ---------------------------------------------------------------------- */
gall_t *gall_build(const tsk_treeseq_t *ts,
                   const int *hap_to_indiv, int n, int use_span)
{
    int S_orig = (int)tsk_treeseq_get_num_trees(ts);

    gall_t *g   = calloc(1, sizeof(*g));
    double *G   = calloc((size_t)n * n, sizeof(double));
    double *in  = calloc((size_t)n, sizeof(double));
    double *out = calloc((size_t)n, sizeof(double));

    tree_cov_matvec_inc_state_t inc;
    int inc_ok = (tree_cov_matvec_inc_state_init(&inc, ts, n, hap_to_indiv) == 0);

    if (!g || !G || !in || !out || !inc_ok) {
        fprintf(stderr, "gall_build: out of memory (n=%d)\n", n);
        free(g); free(G); free(in); free(out);
        if (inc_ok) tree_cov_matvec_inc_state_free(&inc);
        return NULL;
    }

    g->n = n;

    /* Compute the span-weighted eGRM G_all = Σ_s span_s G_s one column at a time
     * using unit-vector inputs (gall_acc_visit applies the per-tree span weight).
     * Cost: n tree-sequence passes, each O(S_orig × n_hap × depth).
     * One-time build; per-phenotype amortises to near-zero. */
    fprintf(stderr, "gall_build: %d tree passes (n=%d, S=%d) ...\n",
            n, n, S_orig);

    gall_acc_ctx_t ctx = { .G = G, .n = n };
    for (int j = 0; j < n; j++) {
        if (j > 0) in[j-1] = 0.0;
        in[j] = 1.0;
        ctx.col = j;
        int ret = tree_walker_run_inc(&inc, S_orig, /*prune_map=*/NULL,
                                      in, out, NULL, gall_acc_visit, &ctx, use_span);
        if (ret != 0) {
            fprintf(stderr, "gall_build: tree walk failed at column %d\n", j);
            goto fail;
        }
    }
    tree_cov_matvec_inc_state_free(&inc);
    free(in);  in  = NULL;
    free(out); out = NULL;

    /* Symmetrise: G_s is symmetric by construction; enforce exact symmetry
     * to guarantee dsyev receives a valid symmetric input. */
    for (int j = 0; j < n; j++) {
        for (int i = j + 1; i < n; i++) {
            double avg = 0.5 * (G[i + (size_t)j*n] + G[j + (size_t)i*n]);
            G[i + (size_t)j*n] = G[j + (size_t)i*n] = avg;
        }
    }

    /* Allocate output arrays and pre-allocated workspace for precond_apply. */
    g->U      = malloc((size_t)n * n * sizeof(double));
    g->lambda = malloc((size_t)n * sizeof(double));
    g->w_ones = malloc((size_t)n * sizeof(double));
    g->work   = malloc((size_t)n * sizeof(double));   /* scratch for precond_apply */
    if (!g->U || !g->lambda || !g->w_ones || !g->work) {
        fprintf(stderr, "gall_build: out of memory for eigen arrays\n");
        goto fail;
    }

    /* Copy G → U (dsyev overwrites input with eigenvectors). */
    memcpy(g->U, G, (size_t)n * n * sizeof(double));
    free(G); G = NULL;

    /* Query workspace, then eigendecompose. */
    fprintf(stderr, "gall_build: eigendecomposing %d×%d G_all ...\n", n, n);
    {
        int lwork = -1, info = 0;
        double lwork_query;
        dsyev_("V", "L", &n, g->U, &n, g->lambda, &lwork_query, &lwork, &info);
        lwork = (int)lwork_query;

        double *dsyev_work = malloc((size_t)lwork * sizeof(double));
        if (!dsyev_work) {
            fprintf(stderr, "gall_build: out of memory for dsyev workspace\n");
            goto fail;
        }
        dsyev_("V", "L", &n, g->U, &n, g->lambda, dsyev_work, &lwork, &info);
        free(dsyev_work);

        if (info != 0) {
            fprintf(stderr, "gall_build: dsyev failed (info=%d)\n", info);
            goto fail;
        }
    }

    /* G_all is PSD by construction (a nonnegative sum of rank-1 branch outer
     * products), so any negative eigenvalue is numerical, not structural.
     * Floor to a small positive value relative to lambda_max so downstream
     * consumers that assume lambda >= 0 stay well-defined:
     * emma_profile_reml() evaluates log(delta*lambda + 1), and
     * gall_precond_build() forms 1/(ws2_mean*lambda + sigma2).  Report when
     * flooring is needed -- a *large* negative lambda signals upstream
     * corruption (see issue #26: OpenMP readout race), not rounding. */
    {
        double lam_max     = g->lambda[n-1];              /* dsyev: ascending */
        double lam_min_raw = g->lambda[0];
        double floor       = 1e-9 * (lam_max > 0.0 ? lam_max : 1.0);
        int    n_floored   = 0;
        for (int i = 0; i < n; i++)
            if (g->lambda[i] < floor) { g->lambda[i] = floor; n_floored++; }
        if (n_floored > 0)
            fprintf(stderr, "gall_build: floored %d non-PSD eigenvalue(s) "
                    "(raw lambda_min = %.3g, floor = %.3g)\n",
                    n_floored, lam_min_raw, floor);
    }

    /* w_ones[k] = (U^T 1_n)[k] = sum of k-th eigenvector column. */
    for (int k = 0; k < n; k++) {
        const double *Uk = g->U + (size_t)k * n;
        double s = 0.0;
        for (int i = 0; i < n; i++) s += Uk[i];
        g->w_ones[k] = s;
    }

    fprintf(stderr, "gall_build: done (λ ∈ [%.3g, %.3g])\n",
            g->lambda[0], g->lambda[n-1]);
    return g;

fail:
    free(G); free(in); free(out);
    if (inc_ok) tree_cov_matvec_inc_state_free(&inc);
    gall_free(g);
    return NULL;
}

void gall_free(gall_t *g)
{
    if (!g) return;
    free(g->U);
    free(g->lambda);
    free(g->w_ones);
    free(g->Pinv);
    free(g->work);
    free(g);
}

/* -------------------------------------------------------------------------
 * File I/O.
 * ---------------------------------------------------------------------- */
int gall_save(const gall_t *g, const char *path, long trees_size)
{
    FILE *f = fopen(path, "wb");
    if (!f) {
        fprintf(stderr, "gall_save: cannot open %s\n", path);
        return -1;
    }

    uint64_t magic = GALL_MAGIC;
    uint32_t ver   = GALL_VERSION;
    int32_t  n32   = g->n;
    int64_t  ts_sz = (int64_t)trees_size;

    int ok = 1;
    ok &= (fwrite(&magic, 8, 1, f) == 1);
    ok &= (fwrite(&ver,   4, 1, f) == 1);
    ok &= (fwrite(&n32,   4, 1, f) == 1);
    ok &= (fwrite(&ts_sz, 8, 1, f) == 1);
    ok &= (fwrite(g->U,      sizeof(double), (size_t)g->n * g->n, f)
               == (size_t)g->n * g->n);
    ok &= (fwrite(g->lambda, sizeof(double), (size_t)g->n, f)
               == (size_t)g->n);
    ok &= (fwrite(g->w_ones, sizeof(double), (size_t)g->n, f)
               == (size_t)g->n);

    fclose(f);
    if (!ok) { fprintf(stderr, "gall_save: write error on %s\n", path); return -1; }
    return 0;
}

gall_t *gall_load(const char *path, int n_expected, long trees_size_expected)
{
    FILE *f = fopen(path, "rb");
    if (!f) return NULL;

    uint64_t magic; uint32_t ver; int32_t n32; int64_t ts_sz;
    int ok = 1;
    ok &= (fread(&magic, 8, 1, f) == 1);
    ok &= (fread(&ver,   4, 1, f) == 1);
    ok &= (fread(&n32,   4, 1, f) == 1);
    ok &= (fread(&ts_sz, 8, 1, f) == 1);

    if (!ok || magic != GALL_MAGIC || ver != GALL_VERSION) {
        fprintf(stderr, "gall_load: bad header in %s\n", path);
        fclose(f); return NULL;
    }
    if ((int)n32 != n_expected) {
        fprintf(stderr, "gall_load: n mismatch (file=%d expected=%d)\n",
                (int)n32, n_expected);
        fclose(f); return NULL;
    }
    if ((long)ts_sz != trees_size_expected) {
        fprintf(stderr, "gall_load: .trees file changed; will rebuild .gall\n");
        fclose(f); return NULL;
    }

    int n = (int)n32;
    gall_t *g = calloc(1, sizeof(*g));
    if (!g) { fclose(f); return NULL; }
    g->n      = n;
    g->U      = malloc((size_t)n * n * sizeof(double));
    g->lambda = malloc((size_t)n * sizeof(double));
    g->w_ones = malloc((size_t)n * sizeof(double));
    g->work   = malloc((size_t)n * sizeof(double));
    if (!g->U || !g->lambda || !g->w_ones || !g->work) goto fail;

    ok  = (fread(g->U,      sizeof(double), (size_t)n * n, f) == (size_t)n * n);
    ok &= (fread(g->lambda, sizeof(double), (size_t)n,     f) == (size_t)n);
    ok &= (fread(g->w_ones, sizeof(double), (size_t)n,     f) == (size_t)n);
    fclose(f);

    if (!ok) {
        fprintf(stderr, "gall_load: read error in %s\n", path);
        gall_free(g); return NULL;
    }
    return g;

fail:
    fclose(f); gall_free(g); return NULL;
}

/* -------------------------------------------------------------------------
 * Preconditioner.
 * ---------------------------------------------------------------------- */
void gall_precond_build(gall_t *g, double ws2_mean, double sigma2)
{
    int n = g->n;
    if (!g->Pinv) {
        g->Pinv = malloc((size_t)n * sizeof(double));
        if (!g->Pinv) return;
    }
    for (int i = 0; i < n; i++)
        g->Pinv[i] = 1.0 / (ws2_mean * g->lambda[i] + sigma2);
    g->ws2_prec  = ws2_mean;
    g->sig2_prec = sigma2;
}

/* P⁻¹ in = U diag(Pinv) Uᵀ in.
 *
 * Pass 1: t[k] = Pinv[k] × (Uᵀ in)[k]    (stored in g->work)
 * Pass 2: out[i] = (U t)[i]
 *
 * g->work is pre-allocated at build/load time; no per-call malloc needed. */
void gall_precond_apply(const gall_t *g, double *out, const double *in)
{
    int n = g->n;
    const double *U    = g->U;
    const double *Pinv = g->Pinv;
    double       *t    = g->work;

    /* t = Pinv × (Uᵀ in) */
    for (int k = 0; k < n; k++) {
        const double *Uk = U + (size_t)k * n;
        double s = 0.0;
        for (int i = 0; i < n; i++) s += Uk[i] * in[i];
        t[k] = Pinv[k] * s;
    }

    /* out = U t */
    for (int i = 0; i < n; i++) {
        double s = 0.0;
        for (int k = 0; k < n; k++) s += U[i + (size_t)k*n] * t[k];
        out[i] = s;
    }
}

/* -------------------------------------------------------------------------
 * CG-compatible wrapper.
 * ---------------------------------------------------------------------- */
void gall_precond_cb(double *out, const double *in, int n, void *data)
{
    (void)n;  /* gall already knows its n */
    gall_precond_apply((gall_t *)data, out, in);
}

/* -------------------------------------------------------------------------
 * EMMA REML 1-D profile-likelihood search.
 *
 * Projects out the intercept exactly (REML), which ensures a well-defined
 * interior maximum even when the ML boundary is flat.
 *
 *   d_i = δ λ_i + 1
 *   A = Σ z_i²/d_i,  B = Σ z_i w_i/d_i,  C = Σ w_i²/d_i
 *   ss_reml = A − B²/C
 *   σ²* = ss_reml/(n−1)
 *   LL = −(n−1)/2 log σ²* − ½ Σ log d_i − ½ log C
 * ---------------------------------------------------------------------- */
static double emma_profile_reml(double log_delta, int n,
                                 const double *lambda,
                                 const double *z, const double *w)
{
    double delta = exp(log_delta);
    double A = 0.0, B = 0.0, C = 0.0, ld = 0.0;
    for (int i = 0; i < n; i++) {
        double d = delta * lambda[i] + 1.0;
        A  += z[i] * z[i] / d;
        B  += z[i] * w[i] / d;
        C  += w[i] * w[i] / d;
        ld += log(d);
    }
    double C_safe   = fmax(C, 1e-30);
    double ss_reml  = A - B * B / C_safe;
    double sig2     = fmax(ss_reml / (n - 1), 1e-15);
    return -(double)(n - 1) * 0.5 * log(sig2) - 0.5 * ld - 0.5 * log(C_safe);
}

int gall_emma_reml(const gall_t *g, const double *y,
                   double *ws2_out, double *sigma2_out)
{
    int n = g->n;
    const double *U      = g->U;
    const double *lambda = g->lambda;
    const double *w      = g->w_ones;   /* U^T 1_n, precomputed */

    /* z = Uᵀ y */
    double *z = malloc((size_t)n * sizeof(double));
    if (!z) { *ws2_out = 0.0; *sigma2_out = 1.0; return 0; }

    for (int k = 0; k < n; k++) {
        const double *Uk = U + (size_t)k * n;
        double s = 0.0;
        for (int i = 0; i < n; i++) s += Uk[i] * y[i];
        z[k] = s;
    }

    /* Golden-section search over log δ ∈ [-20, 20].
     * φ = (√5 − 1)/2 ≈ 0.618 is the golden ratio conjugate. */
    const double PHI = 0.6180339887498949;
    double lo = -20.0, hi = 20.0;
    double x1 = lo + PHI * (hi - lo);
    double x2 = hi - PHI * (hi - lo);
    double f1 = emma_profile_reml(x1, n, lambda, z, w);
    double f2 = emma_profile_reml(x2, n, lambda, z, w);

    for (int iter = 0; iter < 100; iter++) {
        if (f1 > f2) {
            /* x1 (right probe) is better → max is in [x2, hi]; discard [lo, x2] */
            lo = x2; x2 = x1; f2 = f1;
            x1 = lo + PHI * (hi - lo);
            f1 = emma_profile_reml(x1, n, lambda, z, w);
        } else {
            /* x2 (left probe) is better → max is in [lo, x1]; discard [x1, hi] */
            hi = x1; x1 = x2; f1 = f2;
            x2 = hi - PHI * (hi - lo);
            f2 = emma_profile_reml(x2, n, lambda, z, w);
        }
    }

    double delta_opt = exp(0.5 * (x1 + x2));

    /* REML σ²* and ws2* at optimal δ */
    double A = 0.0, B = 0.0, C = 0.0;
    for (int i = 0; i < n; i++) {
        double d = delta_opt * lambda[i] + 1.0;
        A += z[i] * z[i] / d;
        B += z[i] * w[i] / d;
        C += w[i] * w[i] / d;
    }
    double C_safe   = fmax(C, 1e-30);
    double ss_reml  = A - B * B / C_safe;
    *sigma2_out = fmax(ss_reml / (n - 1), 1e-15);
    *ws2_out    = delta_opt * (*sigma2_out);

    free(z);
    return 100;
}
