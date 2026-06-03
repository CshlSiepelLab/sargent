/*
 * SARGENT: Sparse ARG-based ENcoding of Traits
 *
 * Copyright (c) 2026, Cold Spring Harbor Laboratory
 * All rights reserved.
 *
 * This file is part of SARGENT and is distributed under the BSD 3-Clause License.
 * See the LICENSE file in the project root for details.
 */

#pragma once

#include <tskit.h>

/*
 * G_all = Σ_s G_s: aggregate kinship across all local trees.
 *
 * Precomputed once per tree sequence (saved as <prefix>.gall).  Enables two
 * key optimisations in Stage 2:
 *
 *   1. EMMA-based Newton-REML for --uniform-arg mode: replaces 100 EM
 *      iterations (each costing ~700 tree-sequence passes) with a 1-D
 *      profile-likelihood golden-section search — O(100 × n) after the
 *      one-time O(n²) precompute.  Reduces Stage-2 from >60 min to <1 s
 *      for n=1000.
 *
 *   2. G_all preconditioner for per-tree CG: P = ws2_mean * G_all + σ²I.
 *      Because G_all = Σ_s G_s and ws2[s] ≈ ws2_mean on average,
 *      κ(P⁻¹V) ≈ O(1), reducing CG iteration count from ~64 to ~5.
 *
 * File format (.gall):
 *   uint64  magic   = 0x4C4C41475247_4153  ('SARGGALL')
 *   uint32  version = 2
 *   int32   n
 *   int64   trees_size  (st_size of .trees file at build time; used for
 *                        cache invalidation — reload if size changes)
 *   double  U[n*n]      (eigenvectors, col-major; U[:,k] = k-th eigenvector)
 *   double  lambda[n]   (eigenvalues of G_all, ascending order)
 *   double  w_ones[n]   (= U^T 1_n; used for REML intercept projection)
 */

#define GALL_MAGIC   UINT64_C(0x4C4C414753524741)
#define GALL_VERSION 2u

typedef struct {
    int     n;          /* number of individuals */
    double *U;          /* n×n eigenvectors, col-major */
    double *lambda;     /* n eigenvalues, ascending */
    double *w_ones;     /* n entries: (U^T 1_n)[k]; used in REML profile LL */
    double *Pinv;       /* n entries: 1/(ws2_mean*lambda[i]+sigma2); NULL until
                           gall_precond_build() */
    double *work;       /* n scratch doubles for gall_precond_apply(); pre-alloc'd
                           at build/load time to avoid per-call malloc */
    double  ws2_prec;   /* ws2_mean used when Pinv was last built */
    double  sig2_prec;  /* sigma2 used when Pinv was last built */
} gall_t;

/*
 * Build G_all by making n passes through the tree sequence (one per individual
 * basis vector), then eigendecompose via LAPACK dsyev.
 *
 * hap_to_indiv: sample-rank → individual mapping (length = 2*n).
 * Returns NULL on OOM or LAPACK failure.  Caller frees with gall_free().
 */
gall_t *gall_build(const tsk_treeseq_t *ts,
                   const int *hap_to_indiv, int n_indiv);

void gall_free(gall_t *g);

/*
 * Save/load .gall binary file.
 * trees_size: st_size of the corresponding .trees file (used for cache
 *             invalidation in gall_load).
 * gall_save() returns 0 on success, -1 on I/O error.
 * gall_load() returns NULL on I/O error, format mismatch (n != n_expected),
 *             or stale cache (trees_size != trees_size_expected).
 */
int    gall_save(const gall_t *g, const char *path, long trees_size);
gall_t *gall_load(const char *path, int n_expected, long trees_size_expected);

/*
 * Build (or rebuild) the diagonal preconditioner:
 *   Pinv[i] = 1 / (ws2_mean * lambda[i] + sigma2)
 * Safe to call repeatedly when ws2_mean or sigma2 change between CG solves.
 */
void gall_precond_build(gall_t *g, double ws2_mean, double sigma2);

/*
 * Apply preconditioner: out = P⁻¹ in = U diag(Pinv) Uᵀ in.
 * out and in must not alias.  Requires prior gall_precond_build() call.
 * Cost: 2 × n² multiply-adds (two matrix-vector products).
 */
void gall_precond_apply(const gall_t *g, double *out, const double *in);

/*
 * CG-compatible wrapper for gall_precond_apply.
 *
 * Use as cg_opts_t.precond_fn; pass gall_t* as cg_opts_t.precond_data.
 * This satisfies the cg_precond_fn signature:
 *   void f(double *out, const double *in, int n, void *data)
 */
void gall_precond_cb(double *out, const double *in, int n, void *data);

/*
 * EMMA REML profile-likelihood search for the uniform ARG model:
 *
 *   V = sigma2 * I + ws2_total * G_all,   G_all = U Λ Uᵀ
 *
 * Optimises the REML log-likelihood (intercept projected out) by a
 * golden-section search over δ = ws2_total / sigma2 ∈ [e⁻²⁰, e²⁰]:
 *
 *   w = U^T 1_n  (precomputed as g->w_ones)
 *   d_i = δ λ_i + 1
 *   A = Σ z_i²/d_i,  B = Σ z_i w_i/d_i,  C = Σ w_i²/d_i
 *   ss_reml = A − B²/C
 *   σ²* = ss_reml/(n−1)
 *   LL_REML = −(n−1)/2 log σ²* − ½ Σ log d_i − ½ log C
 *
 * On return: *ws2_out = δ* × σ²*, *sigma2_out = σ²*.
 * Returns the number of golden-section iterations (always 100).
 */
int gall_emma_reml(const gall_t *g, const double *y,
                   double *ws2_out, double *sigma2_out);
