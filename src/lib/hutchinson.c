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

#include "hutchinson.h"
#include "geno.h"       /* tree_cov_matvec */

#define HUTCH_DEFAULT_NPROBE 50

static void fill_rademacher(double *v, int n, unsigned int *rng_state)
{
    for (int i = 0; i < n; i++) {
        /* xorshift32 for speed */
        *rng_state ^= *rng_state << 13;
        *rng_state ^= *rng_state >> 17;
        *rng_state ^= *rng_state << 5;
        v[i] = (*rng_state & 1u) ? 1.0 : -1.0;
    }
}

static void fill_gaussian(double *v, int n, unsigned int *rng_state)
{
    /* Box-Muller, consuming pairs */
    for (int i = 0; i < n - 1; i += 2) {
        double u1, u2;
        do {
            *rng_state ^= *rng_state << 13;
            *rng_state ^= *rng_state >> 17;
            *rng_state ^= *rng_state << 5;
            u1 = (*rng_state & 0xFFFFu) / 65536.0;
            *rng_state ^= *rng_state << 13;
            *rng_state ^= *rng_state >> 17;
            *rng_state ^= *rng_state << 5;
            u2 = (*rng_state & 0xFFFFu) / 65536.0;
        } while (u1 < 1e-10);
        double r = sqrt(-2.0 * log(u1));
        double a = 6.28318530717958647693 * u2;
        v[i]   = r * cos(a);
        v[i+1] = r * sin(a);
    }
    if (n & 1) {
        /* handle odd n with a single standard normal via a simple method */
        unsigned int s = *rng_state;
        double u1;
        do {
            s ^= s << 13; s ^= s >> 17; s ^= s << 5;
            u1 = (s & 0xFFFFu) / 65536.0;
        } while (u1 < 1e-10);
        s ^= s << 13; s ^= s >> 17; s ^= s << 5;
        double u2 = (s & 0xFFFFu) / 65536.0;
        v[n-1] = sqrt(-2.0 * log(u1)) * cos(6.28318530717958647693 * u2);
        *rng_state = s;
    }
}

/* dot product of two length-n vectors */
static double vdot(const double *a, const double *b, int n)
{
    double s = 0.0;
    for (int i = 0; i < n; i++) s += a[i] * b[i];
    return s;
}

/* One batch: draw batch_size probes, solve V z = eps for each, then walk the
 * tree sequence accumulating B[s] estimates.
 *
 * eps_buf, z_buf: pre-allocated [batch_size * n] matrices (row-major).
 * tree_work:      pre-allocated [2 * n_nodes] scratch for tree_cov_matvec.
 * gv_buf:         pre-allocated [n] scratch for G_s z_i products.
 */
static void run_batch(cg_state_t *cg, cg_matvec_fn Av, void *av_data,
                      const tsk_treeseq_t *ts, int n, int S,
                      int batch_size, int use_gaussian,
                      unsigned int *rng_state,
                      double *eps_buf, double *z_buf,
                      double *tree_work, tsk_id_t *order_buf,
                      double *gv_buf,
                      double *B, double *B0_accum, double scale)
{
    /* Step 1: draw probe vectors */
    for (int k = 0; k < batch_size; k++) {
        double *eps_k = eps_buf + (size_t)k * n;
        if (use_gaussian)
            fill_gaussian(eps_k, n, rng_state);
        else
            fill_rademacher(eps_k, n, rng_state);
    }

    /* Step 2: solve V z_k = eps_k for each probe.
     * Also accumulate B0 = tr(V^{-1}) = mean(eps_k^T z_k) if requested. */
    for (int k = 0; k < batch_size; k++) {
        double *eps_k = eps_buf + (size_t)k * n;
        double *z_k   = z_buf   + (size_t)k * n;
        memset(z_k, 0, (size_t)n * sizeof(double));  /* cold start */
        cg_solve(cg, Av, av_data, eps_k, z_k);
        if (B0_accum)
            *B0_accum += scale * vdot(eps_k, z_k, n);
    }

    /* Step 3: walk tree sequence; for each tree s accumulate B[s] */
    tsk_tree_t tree;
    int ret = tsk_tree_init(&tree, ts, TSK_NO_SAMPLE_COUNTS);
    if (ret < 0) {
        fprintf(stderr, "hutchinson: tsk_tree_init: %s\n", tsk_strerror(ret));
        return;
    }

    int s = 0;
    for (ret = tsk_tree_first(&tree); ret == TSK_TREE_OK;
         ret = tsk_tree_next(&tree), s++) {
        if (s >= S) break;

        for (int k = 0; k < batch_size; k++) {
            const double *eps_k = eps_buf + (size_t)k * n;
            const double *z_k   = z_buf   + (size_t)k * n;

            /* gv = G_s z_k  (O(n) tree traversal) */
            tree_cov_matvec(gv_buf, z_k, n, &tree, tree_work, order_buf);

            /* accumulate: B[s] += eps_k^T (G_s z_k) */
            B[s] += scale * vdot(eps_k, gv_buf, n);
        }
    }
    if (ret < 0)
        fprintf(stderr, "hutchinson: tree iteration error: %s\n", tsk_strerror(ret));

    tsk_tree_free(&tree);
}

void hutch_accumulate(const double *eps, const double *z, int nprobe, int n,
                      const tsk_treeseq_t *ts, double *B, int S,
                      double *B0_out)
{
    tsk_size_t n_nodes  = tsk_treeseq_get_num_nodes(ts);
    double    *tree_work = malloc(2 * n_nodes * sizeof(double));
    tsk_id_t  *order_buf = malloc(n_nodes * sizeof(tsk_id_t));
    double    *gv_buf    = malloc((size_t)n * sizeof(double));
    if (!tree_work || !order_buf || !gv_buf) {
        fprintf(stderr, "hutch_accumulate: out of memory\n");
        free(tree_work); free(order_buf); free(gv_buf);
        return;
    }

    for (int s = 0; s < S; s++) B[s] = 0.0;
    if (B0_out) *B0_out = 0.0;

    /* B0 = tr(V^{-1}) = mean_k(eps_k^T z_k); computed without any tree walk */
    if (B0_out) {
        for (int k = 0; k < nprobe; k++)
            *B0_out += vdot(eps + (size_t)k * n, z + (size_t)k * n, n);
        *B0_out /= nprobe;
    }

    tsk_tree_t tree;
    int ret = tsk_tree_init(&tree, ts, TSK_NO_SAMPLE_COUNTS);
    if (ret < 0) {
        fprintf(stderr, "hutch_accumulate: tsk_tree_init: %s\n", tsk_strerror(ret));
        free(tree_work); free(order_buf); free(gv_buf);
        return;
    }

    int s = 0;
    for (ret = tsk_tree_first(&tree); ret == TSK_TREE_OK;
         ret = tsk_tree_next(&tree), s++) {
        if (s >= S) break;
        for (int k = 0; k < nprobe; k++) {
            const double *eps_k = eps + (size_t)k * n;
            const double *z_k   = z   + (size_t)k * n;
            tree_cov_matvec(gv_buf, z_k, n, &tree, tree_work, order_buf);
            B[s] += vdot(eps_k, gv_buf, n);
        }
        B[s] /= nprobe;
    }
    if (ret < 0)
        fprintf(stderr, "hutch_accumulate: tree error: %s\n", tsk_strerror(ret));

    tsk_tree_free(&tree);
    free(tree_work);
    free(order_buf);
    free(gv_buf);
}

void hutch_trace_all(cg_state_t *cg, cg_matvec_fn Av, void *av_data,
                     const tsk_treeseq_t *ts, int n, int S,
                     double *B, double *B0_out, const hutch_opts_t *opts)
{
    /* resolve options */
    hutch_opts_t o = {
        .nprobe       = HUTCH_DEFAULT_NPROBE,
        .batch_size   = 0,
        .use_gaussian = 0,
        .seed         = 0,
    };
    if (opts) o = *opts;
    if (o.nprobe <= 0) o.nprobe = HUTCH_DEFAULT_NPROBE;

    int batch = (o.batch_size > 0 && o.batch_size < o.nprobe)
                ? o.batch_size : o.nprobe;

    unsigned int rng = (o.seed != 0) ? o.seed : 0x87654321u;
    /* mix in current rand() state if seed is 0 */
    if (o.seed == 0) rng ^= (unsigned int)rand();

    /* allocate probe and solution matrices for one batch */
    double *eps_buf = malloc((size_t)batch * n * sizeof(double));
    double *z_buf   = malloc((size_t)batch * n * sizeof(double));

    tsk_size_t n_nodes  = tsk_treeseq_get_num_nodes(ts);
    double    *tree_work = malloc(2 * n_nodes * sizeof(double));
    tsk_id_t  *order_buf = malloc(n_nodes * sizeof(tsk_id_t));
    double    *gv_buf    = malloc((size_t)n * sizeof(double));

    if (!eps_buf || !z_buf || !tree_work || !order_buf || !gv_buf) {
        fprintf(stderr, "hutch_trace_all: out of memory\n");
        goto cleanup;
    }

    for (int s = 0; s < S; s++) B[s] = 0.0;
    if (B0_out) *B0_out = 0.0;

    int remaining = o.nprobe;
    while (remaining > 0) {
        int cur = (remaining < batch) ? remaining : batch;
        double scale = 1.0 / o.nprobe;  /* normalise over total probes */

        run_batch(cg, Av, av_data, ts, n, S,
                  cur, o.use_gaussian, &rng,
                  eps_buf, z_buf, tree_work, order_buf, gv_buf,
                  B, B0_out, scale);

        remaining -= cur;
    }

cleanup:
    free(eps_buf);
    free(z_buf);
    free(tree_work);
    free(order_buf);
    free(gv_buf);
}
