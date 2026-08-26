/*
 * SARGENT: Sparse ARG-based ENcoding of Traits
 *
 * Copyright (c) 2026, Cold Spring Harbor Laboratory
 * All rights reserved.
 *
 * This file is part of SARGENT and is distributed under the BSD 3-Clause License.
 * See the LICENSE file in the project root for details.
 *
 * branchsel.c — marginal LMM score over ARG branches, and top-k selection.
 * See include/branchsel.h for the statistic and for why this runs on the full
 * (train + test) ARG rather than the train-only one.
 */

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <stdint.h>
#include <time.h>

#include <tskit.h>
#include "branchsel.h"

/* Portable BLAS declaration (Fortran ABI; same convention as gall.c/predict.c).
 * dsymm: C := alpha * A * B + beta * C with A symmetric. */
extern void dsymm_(const char *side, const char *uplo, const int *m,
                   const int *n, const double *alpha, const double *a,
                   const int *lda, const double *b, const int *ldb,
                   const double *beta, double *c, const int *ldc);

/* Columns of x_b pushed through P per BLAS call.  Bounded so the scratch stays
 * in the tens of MB at any n; the flop count is unaffected. */
#define BS_BLOCK 256

/* ------------------------------------------------------------------------- *
 * Carrier-set dedup.
 *
 * The same clade recurs in every tree it spans, so the raw (tree, node) list is
 * an order of magnitude larger than the set of distinct covariates.  The Python
 * prototype dedups on the exact carrier vector; we do the same, with FNV-1a as
 * the hash and a byte compare to settle collisions.
 * ------------------------------------------------------------------------- */

uint64_t branchsel_fingerprint(const unsigned char *carriers, int n)
{
    uint64_t h = 14695981039346656037ULL;         /* FNV-1a 64 offset basis */
    for (int i = 0; i < n; i++) {
        h ^= (uint64_t)carriers[i];
        h *= 1099511628211ULL;                    /* FNV prime */
    }
    return h;
}

typedef struct {
    uint64_t      *hash;      /* [cap] 0 = empty slot */
    int           *idx;       /* [cap] index into the column store */
    size_t         cap;       /* power of two */
    size_t         used;

    unsigned char *cols;      /* [n_cols * stride] distinct carrier vectors */
    uint64_t      *col_hash;  /* [n_cols] */
    size_t         n_cols;
    size_t         cols_cap;
    int            stride;    /* = n_total */
} bs_table_t;

static int bs_table_init(bs_table_t *t, int stride, size_t hint)
{
    memset(t, 0, sizeof(*t));
    t->stride = stride;
    t->cap = 1024;
    while (t->cap < hint * 2) t->cap <<= 1;
    t->hash = calloc(t->cap, sizeof(uint64_t));
    t->idx  = malloc(t->cap * sizeof(int));
    t->cols_cap = 4096;
    t->cols     = malloc(t->cols_cap * (size_t)stride);
    t->col_hash = malloc(t->cols_cap * sizeof(uint64_t));
    if (!t->hash || !t->idx || !t->cols || !t->col_hash) return -1;
    return 0;
}

static void bs_table_free(bs_table_t *t)
{
    free(t->hash); free(t->idx); free(t->cols); free(t->col_hash);
    memset(t, 0, sizeof(*t));
}

static int bs_table_grow(bs_table_t *t)
{
    size_t ncap = t->cap << 1;
    uint64_t *nh = calloc(ncap, sizeof(uint64_t));
    int      *ni = malloc(ncap * sizeof(int));
    if (!nh || !ni) { free(nh); free(ni); return -1; }
    for (size_t s = 0; s < t->cap; s++) {
        if (t->hash[s] == 0) continue;
        size_t p = t->hash[s] & (ncap - 1);
        while (nh[p] != 0) p = (p + 1) & (ncap - 1);
        nh[p] = t->hash[s];
        ni[p] = t->idx[s];
    }
    free(t->hash); free(t->idx);
    t->hash = nh; t->idx = ni; t->cap = ncap;
    return 0;
}

/* Insert `col` if new.  Returns 0 on success, -1 on OOM. */
static int bs_table_insert(bs_table_t *t, const unsigned char *col)
{
    uint64_t h = branchsel_fingerprint(col, t->stride);
    if (h == 0) h = 1;                            /* 0 marks an empty slot */

    size_t p = h & (t->cap - 1);
    while (t->hash[p] != 0) {
        if (t->hash[p] == h &&
            memcmp(t->cols + (size_t)t->idx[p] * t->stride, col,
                   (size_t)t->stride) == 0)
            return 0;                             /* already have it */
        p = (p + 1) & (t->cap - 1);
    }

    if (t->n_cols == t->cols_cap) {
        size_t nc = t->cols_cap * 2;
        unsigned char *ncols = realloc(t->cols, nc * (size_t)t->stride);
        if (!ncols) return -1;
        t->cols = ncols;
        uint64_t *nhash = realloc(t->col_hash, nc * sizeof(uint64_t));
        if (!nhash) return -1;
        t->col_hash = nhash;
        t->cols_cap = nc;
    }
    memcpy(t->cols + t->n_cols * (size_t)t->stride, col, (size_t)t->stride);
    t->col_hash[t->n_cols] = h;
    t->hash[p] = h;
    t->idx[p]  = (int)t->n_cols;
    t->n_cols++;
    t->used++;

    if (t->used * 2 >= t->cap && bs_table_grow(t) != 0) return -1;
    return 0;
}

/* ------------------------------------------------------------------------- *
 * Tree scan: enumerate branches and collect their distinct carrier vectors.
 * ------------------------------------------------------------------------- */

static int bs_collect(const tsk_treeseq_t *ts, const int *hap_to_indiv,
                      int n_total, int max_trees, bs_table_t *tab,
                      int *n_trees_scanned, int silent)
{
    tsk_size_t nnd    = tsk_treeseq_get_num_nodes(ts);
    tsk_size_t n_hap  = tsk_treeseq_get_num_samples(ts);
    int        S      = (int)tsk_treeseq_get_num_trees(ts);
    const double   *time             = ts->tables->nodes.time;
    const tsk_id_t *sample_index_map = tsk_treeseq_get_sample_index_map(ts);

    /* A tree over n_hap leaves holds at most 2*n_hap - 1 nodes; the slack
     * absorbs multi-root forests and the virtual root if it is ever visited. */
    int maxcol = (int)(2 * n_hap) + 64;

    unsigned char *colbuf   = malloc((size_t)maxcol * n_total);
    int           *col_of   = malloc(nnd * sizeof(int));
    tsk_id_t      *order    = malloc((nnd + 2) * sizeof(tsk_id_t));
    /* keep[s] marks the evenly-spaced subsample; NULL means every tree. */
    unsigned char *keep     = NULL;

    if (!colbuf || !col_of || !order) {
        fprintf(stderr, "branchsel: out of memory (tree scan)\n");
        goto fail;
    }

    if (max_trees > 0 && max_trees < S) {
        keep = calloc((size_t)S, 1);
        if (!keep) goto fail;
        /* linspace(0, S-1, max_trees), rounded — the prototype's subsample. */
        for (int j = 0; j < max_trees; j++) {
            double x = (max_trees == 1) ? 0.0
                     : (double)j * (double)(S - 1) / (double)(max_trees - 1);
            int    t = (int)floor(x + 0.5);
            if (t < 0) t = 0;
            if (t > S - 1) t = S - 1;
            keep[t] = 1;
        }
    }

    tsk_tree_t tree;
    int ret = tsk_tree_init(&tree, ts, 0);
    if (ret < 0) {
        fprintf(stderr, "branchsel: tsk_tree_init: %s\n", tsk_strerror(ret));
        goto fail;
    }

    int scanned = 0;
    int s = 0;
    for (ret = tsk_tree_first(&tree); ret == TSK_TREE_OK;
         ret = tsk_tree_next(&tree), s++) {
        if (s >= S) break;
        if (keep && !keep[s]) continue;
        scanned++;

        for (tsk_size_t i = 0; i < nnd; i++) col_of[i] = -1;

        tsk_size_t n_order = 0;
        if (tsk_tree_postorder(&tree, order, &n_order) < 0) {
            fprintf(stderr, "branchsel: postorder failed at tree %d\n", s);
            tsk_tree_free(&tree);
            goto fail;
        }

        int ncol = 0;
        for (tsk_size_t q = 0; q < n_order; q++) {
            tsk_id_t u = order[q];
            if (u < 0 || (tsk_size_t)u >= nnd) continue;   /* virtual root */
            if (ncol >= maxcol) {
                fprintf(stderr, "branchsel: tree %d exceeded the %d-column "
                        "bound; refusing to truncate\n", s, maxcol);
                tsk_tree_free(&tree);
                goto fail;
            }
            unsigned char *cu = colbuf + (size_t)ncol * n_total;

            if (tree.left_child[u] == TSK_NULL) {
                tsk_id_t si = sample_index_map[u];
                if (si < 0) continue;                      /* non-sample leaf */
                int ind = hap_to_indiv[(int)si];
                if (ind < 0 || ind >= n_total) continue;
                memset(cu, 0, (size_t)n_total);
                cu[ind] = 1;
            } else {
                int first = 1;
                for (tsk_id_t c = tree.left_child[u]; c != TSK_NULL;
                     c = tree.right_sib[c]) {
                    int cc = col_of[c];
                    if (cc < 0) continue;
                    const unsigned char *src = colbuf + (size_t)cc * n_total;
                    if (first) {
                        memcpy(cu, src, (size_t)n_total);
                        first = 0;
                    } else {
                        for (int i = 0; i < n_total; i++) cu[i] += src[i];
                    }
                }
                if (first) continue;                       /* no usable children */
            }
            col_of[u] = ncol;
            ncol++;

            tsk_id_t p = tree.parent[u];
            if (p == TSK_NULL) continue;                   /* root: no branch */
            if (!(time[p] - time[u] > 0.0)) continue;      /* zero-length branch */

            long nc = 0;
            for (int i = 0; i < n_total; i++) nc += cu[i];
            if (nc < 1 || nc > 2L * n_total - 1) continue; /* fixed: no contrast */

            if (bs_table_insert(tab, cu) != 0) {
                fprintf(stderr, "branchsel: out of memory (dedup table)\n");
                tsk_tree_free(&tree);
                goto fail;
            }
        }
    }
    if (ret < 0) {
        fprintf(stderr, "branchsel: tree iteration: %s\n", tsk_strerror(ret));
        tsk_tree_free(&tree);
        goto fail;
    }
    tsk_tree_free(&tree);

    if (!silent)
        fprintf(stderr, "branchsel: %d of %d trees scanned, %zu distinct "
                "branches\n", scanned, S, tab->n_cols);
    *n_trees_scanned = scanned;

    free(colbuf); free(col_of); free(order); free(keep);
    return 0;

fail:
    free(colbuf); free(col_of); free(order); free(keep);
    return -1;
}

/* ------------------------------------------------------------------------- *
 * Scoring and selection.
 * ------------------------------------------------------------------------- */

typedef struct { double z; int j; } bs_rank_t;

static int bs_rank_cmp(const void *a, const void *b)
{
    const bs_rank_t *x = a, *y = b;
    if (x->z > y->z) return -1;
    if (x->z < y->z) return  1;
    return (x->j < y->j) ? -1 : (x->j > y->j);   /* stable, deterministic */
}

int branchsel_run(const tsk_treeseq_t *ts, const int *hap_to_indiv,
                  int n_total, int n_train,
                  const double *P, const double *y_centered,
                  const branchsel_opts_t *opts, branchsel_t *out)
{
    branchsel_opts_t o = { .k = BRANCHSEL_DEFAULT_K,
                           .max_trees = BRANCHSEL_DEFAULT_TREES,
                           .silent = 0 };
    if (opts) o = *opts;
    /* k == 0 is meaningful and must be honoured: it is the no-selection control
     * that has to fall back onto SARGENT-UR exactly.  Only k < 0 means "default". */
    if (o.k < 0) o.k = BRANCHSEL_DEFAULT_K;
    if (o.max_trees < 0) o.max_trees = 0;

    memset(out, 0, sizeof(*out));
    out->n_total = n_total;
    out->n_train = n_train;

    clock_t t0 = clock();

    bs_table_t tab;
    if (bs_table_init(&tab, n_total, 65536) != 0) {
        bs_table_free(&tab);
        return -1;
    }
    int n_trees_scanned = 0;
    if (bs_collect(ts, hap_to_indiv, n_total, o.max_trees, &tab,
                   &n_trees_scanned, o.silent) != 0) {
        bs_table_free(&tab);
        return -1;
    }
    out->n_trees = n_trees_scanned;

    size_t m = tab.n_cols;
    if (m == 0) {
        fprintf(stderr, "branchsel: no usable branches\n");
        bs_table_free(&tab);
        return -1;
    }

    /* Py = P y_c.  P 1 = 0 exactly, so sum(Py) = 0 and the branch mean drops
     * out of the numerator — the raw carrier counts can be used directly. */
    double *Py = malloc((size_t)n_train * sizeof(double));
    double *mu = malloc(m * sizeof(double));
    double *sd = malloc(m * sizeof(double));
    double *num = malloc(m * sizeof(double));
    double *den = malloc(m * sizeof(double));
    double *Cb  = malloc((size_t)n_train * BS_BLOCK * sizeof(double));
    double *PC  = malloc((size_t)n_train * BS_BLOCK * sizeof(double));
    bs_rank_t *rank = malloc(m * sizeof(bs_rank_t));
    if (!Py || !mu || !sd || !num || !den || !Cb || !PC || !rank) {
        fprintf(stderr, "branchsel: out of memory (scoring)\n");
        goto fail_score;
    }

    for (int i = 0; i < n_train; i++) {
        const double *Pi = P + (size_t)i * n_train;
        double acc = 0.0;
        for (int j = 0; j < n_train; j++) acc += Pi[j] * y_centered[j];
        Py[i] = acc;
    }

    /* Train mean / population sd per branch; degenerate branches are dropped. */
    size_t n_keep = 0;
    int *keep_idx = malloc(m * sizeof(int));
    if (!keep_idx) goto fail_score;
    for (size_t j = 0; j < m; j++) {
        const unsigned char *c = tab.cols + j * (size_t)n_total;
        double s1 = 0.0, s2 = 0.0;
        for (int i = 0; i < n_train; i++) { s1 += c[i]; s2 += (double)c[i] * c[i]; }
        double mean = s1 / n_train;
        double var  = s2 / n_train - mean * mean;
        if (var < 0.0) var = 0.0;
        double s = sqrt(var);
        mu[j] = mean; sd[j] = s;
        if (s > 1e-9) keep_idx[n_keep++] = (int)j;
    }
    if (n_keep == 0) {
        fprintf(stderr, "branchsel: every branch is constant within train\n");
        free(keep_idx);
        goto fail_score;
    }

    /* Numerator and denominator, BS_BLOCK branches at a time.
     * den_j = x_j' P x_j, computed exactly (one symmetric matmul per block)
     * rather than by a Hutchinson estimate: at predict-time n the dense P is
     * already formed, and an exact denominator removes probe noise from the
     * ranking entirely. */
    const double one = 1.0, zero = 0.0;
    for (size_t b0 = 0; b0 < n_keep; b0 += BS_BLOCK) {
        int nb = (int)((n_keep - b0 < BS_BLOCK) ? (n_keep - b0) : BS_BLOCK);
        for (int t = 0; t < nb; t++) {
            const unsigned char *c = tab.cols + (size_t)keep_idx[b0 + t] * (size_t)n_total;
            double *dst = Cb + (size_t)t * n_train;    /* column-major for BLAS */
            for (int i = 0; i < n_train; i++) dst[i] = (double)c[i];
        }
        dsymm_("L", "L", &n_train, &nb, &one, P, &n_train,
               Cb, &n_train, &zero, PC, &n_train);
        for (int t = 0; t < nb; t++) {
            const double *cc = Cb + (size_t)t * n_train;
            const double *pc = PC + (size_t)t * n_train;
            double dnum = 0.0, dden = 0.0;
            for (int i = 0; i < n_train; i++) {
                dnum += cc[i] * Py[i];
                dden += cc[i] * pc[i];
            }
            int j = keep_idx[b0 + t];
            num[j] = dnum;
            den[j] = dden;
        }
    }

    for (size_t q = 0; q < n_keep; q++) {
        int j = keep_idx[q];
        double d = den[j];
        rank[q].z = (d > 1e-12) ? fabs(num[j]) / sqrt(d) : 0.0;
        rank[q].j = j;
    }
    qsort(rank, n_keep, sizeof(bs_rank_t), bs_rank_cmp);

    /* Optional per-branch dump, keyed by carrier fingerprint so an external
     * reference implementation can align its branches with ours. */
    {
        const char *dump = getenv("SARGENT_BRANCHSEL_DUMP");
        if (dump && *dump) {
            FILE *df = fopen(dump, "w");
            if (df) {
                fprintf(df, "fingerprint\tn_carrier_total\tn_carrier_train\t"
                            "num\tden\tz\trank\n");
                for (size_t q = 0; q < n_keep; q++) {
                    int j = rank[q].j;
                    const unsigned char *c = tab.cols + (size_t)j * (size_t)n_total;
                    long nct = 0, ncr = 0;
                    for (int i = 0; i < n_total; i++) {
                        nct += c[i];
                        if (i < n_train) ncr += c[i];
                    }
                    fprintf(df, "%016llx\t%ld\t%ld\t%.17g\t%.17g\t%.17g\t%zu\n",
                            (unsigned long long)tab.col_hash[j], nct, ncr,
                            num[j], den[j], rank[q].z, q);
                }
                fclose(df);
                if (!o.silent)
                    fprintf(stderr, "branchsel: wrote per-branch scores to %s\n", dump);
            } else {
                fprintf(stderr, "branchsel: cannot write %s\n", dump);
            }
        }
    }

    int k = (o.k < (int)n_keep) ? o.k : (int)n_keep;
    out->n_branch = (int)n_keep;
    out->k = k;
    out->Z     = (k > 0) ? malloc((size_t)n_total * k * sizeof(double)) : NULL;
    out->score = (k > 0) ? malloc((size_t)k * sizeof(double)) : NULL;
    if (k > 0 && (!out->Z || !out->score)) { free(keep_idx); goto fail_score; }

    for (int t = 0; t < k; t++) {
        int j = rank[t].j;
        const unsigned char *c = tab.cols + (size_t)j * (size_t)n_total;
        double m_j = mu[j], s_j = sd[j];
        for (int i = 0; i < n_total; i++)
            out->Z[(size_t)i * k + t] = ((double)c[i] - m_j) / s_j;
        out->score[t] = rank[t].z;
    }

    out->seconds = (double)(clock() - t0) / CLOCKS_PER_SEC;
    if (!o.silent) {
        if (k > 0)
            fprintf(stderr, "branchsel: scored %d branches, selected top %d "
                    "(z: %.3f .. %.3f) in %.1f s CPU\n",
                    out->n_branch, k, out->score[0], out->score[k - 1],
                    out->seconds);
        else
            fprintf(stderr, "branchsel: scored %d branches, selected none "
                    "(k=0 control) in %.1f s CPU\n",
                    out->n_branch, out->seconds);
    }

    free(keep_idx);
    free(Py); free(mu); free(sd); free(num); free(den);
    free(Cb); free(PC); free(rank);
    bs_table_free(&tab);
    return 0;

fail_score:
    free(Py); free(mu); free(sd); free(num); free(den);
    free(Cb); free(PC); free(rank);
    free(out->Z); free(out->score);
    out->Z = NULL; out->score = NULL;
    bs_table_free(&tab);
    return -1;
}

void branchsel_free(branchsel_t *r)
{
    if (!r) return;
    free(r->Z); free(r->score);
    r->Z = NULL; r->score = NULL;
}

void branchsel_kinship(const branchsel_t *r, double *K_train, double *K_cross)
{
    int n_tr = r->n_train, n_te = r->n_total - r->n_train, k = r->k;
    if (k <= 0) return;             /* no-selection control: caller drops the component */
    double inv_k = 1.0 / (double)k;

    for (int i = 0; i < n_tr; i++) {
        const double *zi = r->Z + (size_t)i * k;
        for (int j = 0; j <= i; j++) {
            const double *zj = r->Z + (size_t)j * k;
            double s = 0.0;
            for (int t = 0; t < k; t++) s += zi[t] * zj[t];
            s *= inv_k;
            K_train[(size_t)i * n_tr + j] = s;
            K_train[(size_t)j * n_tr + i] = s;
        }
    }
    for (int c = 0; c < n_te; c++) {
        const double *zc = r->Z + (size_t)(n_tr + c) * k;
        for (int j = 0; j < n_tr; j++) {
            const double *zj = r->Z + (size_t)j * k;
            double s = 0.0;
            for (int t = 0; t < k; t++) s += zc[t] * zj[t];
            K_cross[(size_t)c * n_tr + j] = s * inv_k;
        }
    }
}
