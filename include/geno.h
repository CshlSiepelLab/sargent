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

#include <stddef.h>
#include <tskit.h>

/* Flags for read_variant / read_chunk */
#define GENO_RAW         0
#define GENO_IMPUTE      (1<<0)   /* replace NaN missing with per-SNP mean */
#define GENO_STANDARDIZE (1<<1)   /* center + scale to variance 1; implies IMPUTE */

typedef struct {
    double mean;    /* E[g] over non-missing calls  (= 2 * alt_freq) */
    double scale;   /* 1/sqrt(2*p*(1-p)); monomorphic variants are excluded */
    int    n_miss;  /* number of missing calls */
} geno_var_stats_t;

typedef struct geno_reader geno_reader_t;

typedef struct {
    /* One streaming pass computing mean/scale for every polymorphic variant.
       Sets r->n_variants to the polymorphic count. Must be called before
       read_variant / read_chunk / seek. */
    int  (*precompute_stats)(geno_reader_t *r);

    /* Decode the next polymorphic variant into out[r->n_samples].
       Returns 1 on success, 0 at EOF, -1 on error. */
    int  (*read_variant)(geno_reader_t *r, double *out, int flags);

    /* Decode up to k polymorphic variants into out[k * r->n_samples] row-major.
       Returns the number of variants decoded (0 = EOF, -1 = error). */
    int  (*read_chunk)(geno_reader_t *r, int k, double *out, int flags);

    /* Seek to the i-th polymorphic variant (0-based). */
    int  (*seek)(geno_reader_t *r, int var_idx);

    void (*rewind)(geno_reader_t *r);
    void (*close)(geno_reader_t *r);
} geno_ops_t;

struct geno_reader {
    const geno_ops_t  *ops;
    int                n_samples;   /* number of individuals */
    int                n_variants;  /* polymorphic variants (set by precompute_stats) */
    geno_var_stats_t  *stats;       /* [n_variants], NULL until precompute_stats */
    void              *impl;        /* backend-private state */
};

#define geno_precompute_stats(r)          (r)->ops->precompute_stats(r)
#define geno_read_variant(r, out, flags)  (r)->ops->read_variant(r, out, flags)
#define geno_read_chunk(r, k, out, flags) (r)->ops->read_chunk(r, k, out, flags)
#define geno_seek(r, i)                   (r)->ops->seek(r, i)
#define geno_rewind(r)                    (r)->ops->rewind(r)
#define geno_close(r)                     (r)->ops->close(r)

/* Backend constructors */
geno_reader_t *geno_open_bed(const char *prefix);  /* .bed/.bim/.fam prefix */

/* Caller owns ts and is responsible for tsk_treeseq_free — the reader borrows it.
   This allows the same tsk_treeseq_t to be used concurrently for local-tree
   covariance matrices and other tree-sequence operations. */
geno_reader_t *geno_open_trees(tsk_treeseq_t *ts);

/* Build a sample-rank → individual-index mapping from the tree sequence's
 * individual table.  Returns a malloc'd array of length n_haplotypes where
 * result[rank] = individual index for haplotype sample rank `rank`.
 *
 * Requires every individual to be diploid (exactly 2 sample nodes).  Falls
 * back to rank/2 with a warning if the tree sequence has no individual table.
 * Returns NULL on allocation failure or inconsistent ploidy.  Caller must free. */
int *tsk_build_hap_to_indiv(const tsk_treeseq_t *ts);

/* Compute out = G_s * in for the local tree pointed to by `tree`, where
 * G_s[i,j] = total branch length shared between diploid individuals i and j
 * (sum of branch lengths from root down to their MRCA).
 *
 * n_indiv:      number of diploid individuals.
 * hap_to_indiv: sample-rank → individual-index mapping, length n_haplotypes.
 *               Build once with tsk_build_hap_to_indiv and reuse across trees.
 * work:         caller-allocated scratch, 2 * tsk_treeseq_get_num_nodes(ts) doubles.
 * order_buf:    caller-allocated scratch, tsk_treeseq_get_num_nodes(ts) tsk_id_t values.
 *
 * work and order_buf should be allocated once and reused across tree iterations. */
void tree_cov_matvec(double *out, const double *in, int n_indiv,
                     const tsk_tree_t *tree,
                     double *work, tsk_id_t *order_buf,
                     const int *hap_to_indiv);
