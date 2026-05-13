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
