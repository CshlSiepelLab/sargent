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
#include <string.h>
#include <math.h>
#include <float.h>

#include "conjgrad.h"

#define CG_DEFAULT_TOL      1e-6
#define CG_DEFAULT_ABS_TOL  0.0

static void opts_defaults(cg_opts_t *o, int n)
{
    if (o->max_iter <= 0) o->max_iter = n;
    if (o->tol <= 0.0)    o->tol      = CG_DEFAULT_TOL;
}

cg_state_t *cg_state_new(int n, const cg_opts_t *opts)
{
    cg_state_t *st = calloc(1, sizeof(*st));
    if (!st) return NULL;

    st->n = n;
    st->r = malloc((size_t)n * sizeof(double));
    st->p = malloc((size_t)n * sizeof(double));
    st->q = malloc((size_t)n * sizeof(double));
    st->z = malloc((size_t)n * sizeof(double));  /* used only with preconditioner */

    if (!st->r || !st->p || !st->q || !st->z) {
        cg_state_free(st);
        return NULL;
    }

    if (opts) {
        st->opts = *opts;
    } else {
        st->opts.tol          = CG_DEFAULT_TOL;
        st->opts.abs_tol      = CG_DEFAULT_ABS_TOL;
        st->opts.max_iter     = 0;   /* resolved to n in opts_defaults */
        st->opts.diag_precond = NULL;
    }
    opts_defaults(&st->opts, n);
    return st;
}

void cg_state_free(cg_state_t *st)
{
    if (!st) return;
    free(st->r);
    free(st->p);
    free(st->q);
    free(st->z);
    free(st);
}

void cg_state_set_opts(cg_state_t *st, const cg_opts_t *opts)
{
    if (!opts) return;
    st->opts = *opts;
    opts_defaults(&st->opts, st->n);
}

/* dot product */
static double dot(const double *a, const double *b, int n)
{
    double s = 0.0;
    for (int i = 0; i < n; i++) s += a[i] * b[i];
    return s;
}

/* y += alpha * x */
static void axpy(double *y, double alpha, const double *x, int n)
{
    for (int i = 0; i < n; i++) y[i] += alpha * x[i];
}

int cg_solve(cg_state_t *st, cg_matvec_fn Av, void *data,
             const double *b, double *x)
{
    int n = st->n;
    double *r = st->r;
    double *p = st->p;
    double *q = st->q;
    double *z = st->z;   /* preconditioned residual */

    const double     *M  = st->opts.diag_precond;   /* NULL → no diagonal precond */
    const cg_precond_fn Pfn  = st->opts.precond_fn;     /* NULL → fall back to diag */
    void             *Pd  = st->opts.precond_data;

    /* r = b - A x */
    Av(q, x, n, data);
    for (int i = 0; i < n; i++) r[i] = b[i] - q[i];

    double b_norm = sqrt(dot(b, b, n));
    if (b_norm < DBL_EPSILON) {
        memset(x, 0, (size_t)n * sizeof(double));
        st->last_iters = 0;
        st->last_resid = 0.0;
        return 0;
    }

    /* z = M⁻¹ r  (general preconditioner takes priority over diagonal) */
    if (Pfn) {
        Pfn(z, r, n, Pd);
    } else if (M) {
        for (int i = 0; i < n; i++) z[i] = M[i] * r[i];
    } else {
        memcpy(z, r, (size_t)n * sizeof(double));
    }

    memcpy(p, z, (size_t)n * sizeof(double));

    double rz = dot(r, z, n);

    int max_iter = st->opts.max_iter;
    double tol   = st->opts.tol;
    double atol  = st->opts.abs_tol;
    int iter;

    for (iter = 0; iter < max_iter; iter++) {
        double r_norm = sqrt(dot(r, r, n));
        double rel    = r_norm / b_norm;

        st->last_resid = rel;
        if (rel < tol) break;
        if (atol > 0.0 && r_norm < atol) break;

        Av(q, p, n, data);

        double pAp = dot(p, q, n);
        if (pAp <= 0.0) break;        /* loss of positive-definiteness */

        double alpha = rz / pAp;

        axpy(x,  alpha,   p, n);
        axpy(r, -alpha,   q, n);

        /* z = M⁻¹ r */
        if (Pfn) {
            Pfn(z, r, n, Pd);
        } else if (M) {
            for (int i = 0; i < n; i++) z[i] = M[i] * r[i];
        } else {
            memcpy(z, r, (size_t)n * sizeof(double));
        }

        double rz_new = dot(r, z, n);
        double beta   = rz_new / rz;

        for (int i = 0; i < n; i++) p[i] = z[i] + beta * p[i];

        rz = rz_new;
    }

    st->last_iters = iter;
    st->last_resid = sqrt(dot(r, r, n)) / b_norm;
    return iter;
}
