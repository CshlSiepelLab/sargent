/*
 * SARGENT: Sparse ARG-based ENcoding of Traits
 *
 * Copyright (c) 2026, Cold Spring Harbor Laboratory
 * All rights reserved.
 *
 * This file is part of SARGENT and is distributed under the BSD 3-Clause License.
 * See the LICENSE file in the project root for details.
 */

/* In-memory genotype reader over a ROW SUBSET of a matrix already held in core.
 *
 * Written for Stage-1 cross-fitting (--stage1-cv).  Cross-fitting has to refit
 * Stage 1 on K-1 folds of the training individuals at a time, and the on-disk
 * readers address whole files, so there was no way to hand var_fixed a subset
 * of rows.  This backend wraps predict.c's in-core X_train and serves only the
 * requested rows, which is enough because Stage 1 touches genotypes solely
 * through geno_read_variant.
 *
 * The matrix handed in is variant-major, X[j * n_full + i], and is already
 * standardized against the FULL training sample.  This reader re-centres and
 * re-scales each variant WITHIN the subset, so a fold's model is fitted in the
 * fold's own units and never sees the moments of the held-out individuals.
 * The resulting per-variant transform is published in r->stats, so the caller
 * can apply exactly the same transform when scoring the held-out rows -- which
 * is what keeps the out-of-fold prediction on the same scale as the fit.
 *
 * A variant that is monomorphic within a subset is NOT dropped: doing so would
 * renumber the variants and break alignment with beta[]. It is served as zeros,
 * so it simply carries no effect for that fold.
 */

#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "geno.h"

typedef struct {
    const double *X;        /* borrowed, variant-major [p * n_full] */
    int           p;
    int           n_full;
    int          *rows;     /* owned copy, [n_rows] indices into 0..n_full-1 */
    int           n_rows;
    int           cur;      /* streaming cursor over variants */
} mem_impl_t;

static int mem_precompute_stats(geno_reader_t *r)
{
    mem_impl_t *m = (mem_impl_t *)r->impl;
    if (r->stats) return 0;                    /* idempotent */

    r->stats = calloc((size_t)m->p, sizeof(geno_var_stats_t));
    if (!r->stats) return -1;

    for (int j = 0; j < m->p; j++) {
        const double *col = m->X + (size_t)j * m->n_full;
        double s = 0.0;
        for (int t = 0; t < m->n_rows; t++) s += col[m->rows[t]];
        double mean = s / m->n_rows;

        double ss = 0.0;
        for (int t = 0; t < m->n_rows; t++) {
            double d = col[m->rows[t]] - mean;
            ss += d * d;
        }
        double var = ss / m->n_rows;

        r->stats[j].mean   = mean;
        /* Monomorphic within this subset: scale 0 serves the variant as zeros
         * rather than renumbering it away, keeping beta[] alignment intact.
         *
         * The test MUST be a tolerance, not var > 0.  Summing (x-mean)^2 over
         * identical values gives ~1e-31 rather than exactly zero, so an exact
         * test lets 1/sqrt(var) reach ~3e15; a held-out carrier of that variant
         * then contributes a term of ~1e15 to its PRS and the out-of-fold
         * vector explodes (measured: var(prs_oof) ~ 1e28 against 0.7 in
         * sample).  The incoming matrix is standardized to unit variance on the
         * FULL training sample, so an in-fold variance this far below 1 means
         * the variant is constant here and carries no information for this
         * fold, whatever the floating-point residue says. */
        r->stats[j].scale  = (var > 1e-8) ? 1.0 / sqrt(var) : 0.0;
        r->stats[j].n_miss = 0;
    }
    r->n_variants = m->p;
    return 0;
}

static int mem_read_at(geno_reader_t *r, int j, double *out, int flags)
{
    mem_impl_t *m = (mem_impl_t *)r->impl;
    if (!r->stats) return -1;
    if (j < 0 || j >= m->p) return -1;

    const double *col = m->X + (size_t)j * m->n_full;
    if (flags & GENO_STANDARDIZE) {
        double mu = r->stats[j].mean, sc = r->stats[j].scale;
        for (int t = 0; t < m->n_rows; t++)
            out[t] = (col[m->rows[t]] - mu) * sc;
    } else {
        for (int t = 0; t < m->n_rows; t++)
            out[t] = col[m->rows[t]];
    }
    return 1;
}

static int mem_read_variant(geno_reader_t *r, double *out, int flags)
{
    mem_impl_t *m = (mem_impl_t *)r->impl;
    if (m->cur >= m->p) return 0;                /* EOF */
    int rc = mem_read_at(r, m->cur, out, flags);
    if (rc > 0) m->cur++;
    return rc;
}

static int mem_read_chunk(geno_reader_t *r, int k, double *out, int flags)
{
    mem_impl_t *m = (mem_impl_t *)r->impl;
    int got = 0;
    for (; got < k && m->cur < m->p; got++) {
        if (mem_read_at(r, m->cur, out + (size_t)got * m->n_rows, flags) <= 0) break;
        m->cur++;
    }
    return got;
}

static int mem_seek(geno_reader_t *r, int var_idx)
{
    mem_impl_t *m = (mem_impl_t *)r->impl;
    if (var_idx < 0 || var_idx >= m->p) return -1;
    m->cur = var_idx;
    return 0;
}

static void mem_rewind(geno_reader_t *r)
{
    ((mem_impl_t *)r->impl)->cur = 0;
}

static void mem_close(geno_reader_t *r)
{
    if (!r) return;
    mem_impl_t *m = (mem_impl_t *)r->impl;
    if (m) { free(m->rows); free(m); }
    free(r->stats);
    free(r);
}

static const geno_ops_t MEM_OPS = {
    .precompute_stats = mem_precompute_stats,
    .read_variant     = mem_read_variant,
    .read_chunk       = mem_read_chunk,
    .seek             = mem_seek,
    .rewind           = mem_rewind,
    .close            = mem_close,
};

geno_reader_t *geno_open_mem_subset(const double *X, int p, int n_full,
                                    const int *rows, int n_rows)
{
    if (!X || !rows || p <= 0 || n_rows <= 1 || n_rows > n_full) return NULL;

    geno_reader_t *r = calloc(1, sizeof(*r));
    mem_impl_t    *m = calloc(1, sizeof(*m));
    if (!r || !m) { free(r); free(m); return NULL; }

    m->rows = malloc((size_t)n_rows * sizeof(int));
    if (!m->rows) { free(r); free(m); return NULL; }
    memcpy(m->rows, rows, (size_t)n_rows * sizeof(int));

    m->X = X; m->p = p; m->n_full = n_full; m->n_rows = n_rows; m->cur = 0;

    r->ops        = &MEM_OPS;
    r->impl       = m;
    r->n_samples  = n_rows;
    r->n_variants = p;      /* provisional; precompute_stats confirms */
    r->stats      = NULL;
    return r;
}
