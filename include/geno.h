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

/* ---------------------------------------------------------------------------
 * BED-only helpers for out-of-sample prediction (native --predict path).
 *
 * A prediction run must align a *test* BED to the *train* polymorphic variant
 * set (β is indexed by that set) and standardize test dosages by TRAIN
 * statistics.  Train and test share the same .bim (same variants, same order),
 * so alignment is by absolute file-variant index.
 * ------------------------------------------------------------------------- */

/* Return the absolute .bim/.bed file-variant index of the j-th polymorphic
 * variant (0 <= j < r->n_variants), i.e. the inverse of the monomorphic
 * filtering done by precompute_stats.  Returns -1 for a non-BED reader or an
 * out-of-range j.  Requires precompute_stats to have been called. */
int geno_bed_poly_file_index(const geno_reader_t *r, int j);

/* Decode the raw dosages (0,1,2,NaN; no imputation/standardization) of the
 * absolute file-variant index `file_var_idx` into out[r->n_samples].  Valid
 * before or after precompute_stats (does not depend on the poly filter).
 * Returns 0 on success, -1 on error / non-BED reader.  Used to read test-BED
 * variants at train-derived file indices. */
int geno_bed_read_raw_at(geno_reader_t *r, int file_var_idx, double *out);

/* Caller owns ts and is responsible for tsk_treeseq_free — the reader borrows it.
   This allows the same tsk_treeseq_t to be used concurrently for local-tree
   covariance matrices and other tree-sequence operations. */
geno_reader_t *geno_open_trees(tsk_treeseq_t *ts);

/* In-memory reader over a ROW SUBSET of an already-loaded, variant-major matrix
 * X[j * n_full + i].  Used by Stage-1 cross-fitting, which must refit on K-1
 * folds of the training individuals and so cannot address whole files.
 *
 * Each variant is re-centred and re-scaled WITHIN the subset, so a fold's model
 * never sees the moments of its held-out individuals; the transform is exposed
 * through r->stats[j].{mean,scale} so the caller can score the held-out rows on
 * the same scale.  A variant that is monomorphic within the subset is served as
 * zeros rather than dropped, which keeps variant indices aligned with beta[].
 * The matrix is BORROWED and must outlive the reader. */
geno_reader_t *geno_open_mem_subset(const double *X, int p, int n_full,
                                    const int *rows, int n_rows);

/* Estimate E[tr(G_s)/n] using a single Rademacher probe over all trees.
 * Used to calibrate the tau2_mix prior grid to phenotype-variance units
 * regardless of branch-length scale (e.g. generations vs coalescent units).
 * Returns 1.0 on allocation failure. */
double estimate_g_scale(const tsk_treeseq_t *ts, const int *hap_to_indiv,
                        int n_indiv, unsigned int rng_seed);

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

/* ---------------------------------------------------------------------------
 * Incremental tree-cov matvec.
 *
 * Maintains S[u] (post-order subtree sum of `in`), our_parent[u], and
 * bl[u] (= time[parent(u)] - time[u]) across tree boundaries using tskit's
 * edge-diff machinery.  Each tree boundary pays O(|edges_changed| * depth)
 * for S[] and our_parent[] / bl[] updates, plus a per-sample walk-up from
 * each sample to the root for the matvec output (no per-tree pre-order
 * pass over n_nodes).
 *
 * Usage:
 *
 *   tree_cov_matvec_inc_state_t st;
 *   tree_cov_matvec_inc_state_init(&st, ts, n_indiv, hap_to_indiv);
 *   for each matvec call:
 *       tree_cov_matvec_inc_state_reset(&st, in);  // resets S[], our_parent[], bl[]
 *       for ret = tsk_tree_first(&tree); ret == OK; ret = tsk_tree_next(&tree):
 *           tree_cov_matvec_incremental(out, &st, &tree);
 *           ...
 *   tree_cov_matvec_inc_state_free(&st);
 *
 * Validated against the naive tree_cov_matvec and the pairwise reference via
 * test/lib/check_tree_matvec.
 * --------------------------------------------------------------------------- */

typedef struct {
    const tsk_treeseq_t *ts;
    int                  n_indiv;
    const int           *hap_to_indiv;   /* not owned */
    tsk_size_t           n_nodes;
    double              *S;              /* length n_nodes; updated incrementally */
    tsk_id_t            *our_parent;     /* length n_nodes; tracks tree structure */
    double              *bl;             /* length n_nodes; bl[u] = time[parent(u)] - time[u] */
    /* Readout scratch.  `time_order` lists every node once, sorted by
     * DECREASING node time.  tskit requires time[parent] > time[child] on
     * every edge, so this single global order is a valid parent-before-child
     * order for *every* tree in the sequence -- no per-tree preorder needed.
     * `cum` holds the root-path prefix sum; see tree_cov_matvec_inc_readout. */
    tsk_id_t            *time_order;     /* length n_nodes; decreasing node time */
    double              *cum;            /* length n_nodes; scratch */

    /* Per-tree compaction, for INFERRED ARGs.
     *
     * `time_order` spans the whole node table, which is the right thing for a
     * simulated ARG: there n_nodes (7,444 on the benchmark cell) barely exceeds
     * the ~5,998 present in any one tree, so sweeping all of them wastes
     * nothing.  A Relate ARG is the opposite shape -- 368,878 nodes over 122
     * trees, still ~5,998 per tree -- so the same sweep did 60x redundant work
     * and ran ~5x SLOWER than the leaf-to-root loop it replaced.
     *
     * When the node table dwarfs a single tree we therefore sweep a compacted,
     * time-ordered list of the nodes actually present.  Building one costs an
     * O(n_nodes) filter, so it is cached per tree index and reused across the
     * n_tr column passes, which all traverse the same sequence. */
    int                  compact;        /* 1 = use the per-tree lists below */
    int                  cur_tree;       /* tree index; set by tree_walker_run_inc */
    tsk_size_t           n_trees;
    tsk_id_t           **tree_nodes;     /* [n_trees], NULL until first visited */
    tsk_size_t          *tree_nnodes;    /* [n_trees] */
    unsigned char       *has_child;      /* n_nodes; scratch for the filter */
} tree_cov_matvec_inc_state_t;

/* Allocate scratch buffers; returns 0 on success, -1 on OOM. */
int  tree_cov_matvec_inc_state_init(tree_cov_matvec_inc_state_t *st,
                                    const tsk_treeseq_t *ts,
                                    int n_indiv,
                                    const int *hap_to_indiv);

void tree_cov_matvec_inc_state_free(tree_cov_matvec_inc_state_t *st);

/* Reset S[] and our_parent[] to the "empty tree" state for the given input
 * vector.  Call once per matvec before iterating trees. */
void tree_cov_matvec_inc_state_reset(tree_cov_matvec_inc_state_t *st,
                                     const double *in);

/* Process the edge diff at tree boundary using tree->tree_pos and compute
 * out = G_s * in for the current tree.  Call AFTER tsk_tree_first/_next has
 * positioned `tree` at the new interval.
 *
 * Internally this is `tree_cov_matvec_inc_advance(st, tree)` followed by
 * `tree_cov_matvec_inc_readout(out, st)`.  Use the split variants directly
 * if you want to skip the readout for trees that don't contribute to the
 * accumulator (e.g. ws²[sp] == 0 in v_matvec). */
void tree_cov_matvec_incremental(double *out,
                                 tree_cov_matvec_inc_state_t *st,
                                 const tsk_tree_t *tree);

/* Apply the edge diff at tree->tree_pos to S[], our_parent[], bl[].  Must be
 * called for every tree boundary even if the corresponding readout is
 * skipped -- otherwise S[] drifts out of sync with the underlying tree. */
void tree_cov_matvec_inc_advance(tree_cov_matvec_inc_state_t *st,
                                  const tsk_tree_t *tree);

/* Per-sample walk-up readout: out[i] = sum over haps h of indiv i of
 * sum_{v ancestor of h} bl[v] * S[v].  Uses the current contents of S[],
 * our_parent[], bl[] in `st` (must be coherent with the current tree). */
void tree_cov_matvec_inc_readout(double *out,
                                  const tree_cov_matvec_inc_state_t *st);

/* Reference implementation of out = G_s * in, computed by direct pairwise
 * MRCA-time summation.
 *
 * For each haploid sample pair (h_a, h_b):
 *   G contribution to (indiv(h_a), indiv(h_b)) =
 *       0.5 * (time(root) - time(MRCA(h_a, h_b)))
 *
 * O(n_haps^2 * depth) per tree.  Intended for the paranoia harness
 * (test/lib/check_tree_matvec) -- cross-check against tree_cov_matvec
 * (current production) and tree_cov_matvec_incremental (Phase 5 chunk 3).
 * Do NOT use in the hot path.
 *
 * mark: caller-allocated scratch, n_nodes ints, contents arbitrary on entry. */
void tree_cov_matvec_pairwise(double *out, const double *in, int n_indiv,
                              const tsk_tree_t *tree,
                              int *mark,
                              const int *hap_to_indiv);

/* ---------------------------------------------------------------------------
 * Tree-walker API (Phase 5).
 *
 * Iterates a tree sequence once, calling a per-tree callback with the result
 * of G_s * in for each tree.  Centralises tree-by-tree iteration across the
 * three callers (v_matvec, compute_quad_all, hutch trace probe).  The
 * implementation can be evolved from naive per-tree traversal to incremental
 * edge-diff updates without touching call sites.
 *
 * Usage:
 *
 *   int my_pre(int s_orig, int sp, void *ctx) { ... return 1 to visit; }
 *   void my_visit(const tree_walk_state_t *st, void *ctx) { ... }
 *   tree_walker_run(ts, hap_to_indiv, n_indiv, S_orig, prune_map,
 *                   in, gv, work, order_buf, my_pre, my_visit, ctx);
 *
 * pre_visit may be NULL (always visit).
 * --------------------------------------------------------------------------- */

typedef struct {
    int      n_indiv;        /* convenience copy */
    int      s_orig;         /* current original-tree index in [0, S_orig) */
    int      sp;             /* current pruned-segment index (= s_orig when prune_map is NULL) */
    double   span;           /* span weight of the current tree = (span_s / mean_span);
                                aggregation visitors weight G_s by this so the kernel
                                is the span-weighted eGRM Sum_s span_s G_s (normalised so
                                Sum_s span == S_orig, preserving the old unweighted scale) */
    double  *gv;             /* G_s * in for the current tree, length n_indiv */
} tree_walk_state_t;

typedef int  (*tree_pre_visit_fn)(int s_orig, int sp, void *ctx);
typedef void (*tree_visit_fn)    (const tree_walk_state_t *st, void *ctx);

/* Returns 0 on success, -1 on tskit init/iteration failure. */
int tree_walker_run(const tsk_treeseq_t *ts,
                    const int           *hap_to_indiv,
                    int                  n_indiv,
                    int                  S_orig,
                    const int           *prune_map,    /* nullable; NULL = identity */
                    const double        *in,
                    double              *gv,            /* length n_indiv */
                    double              *work,          /* length 2 * n_nodes */
                    tsk_id_t            *order_buf,     /* length n_nodes */
                    tree_pre_visit_fn    pre_visit,     /* nullable */
                    tree_visit_fn        visit,
                    void                *ctx);

/* ---------------------------------------------------------------------------
 * Incremental tree-walker (Phase 5 chunk 4): same shape as tree_walker_run
 * but consumes an externally-owned inc_state for the per-tree matvec.
 *
 * - calls tree_cov_matvec_inc_state_reset(inc, in) at the start.
 * - calls tree_cov_matvec_inc_advance(inc, &tree) at every tree boundary
 *   (state must stay coherent even when readout is skipped).
 * - calls pre_visit; if it returns 0, skip the readout but keep state.
 * - otherwise calls tree_cov_matvec_inc_readout(gv, inc) and then visit.
 *
 * Returns 0 on success, -1 on tskit error.
 * --------------------------------------------------------------------------- */
int tree_walker_run_inc(tree_cov_matvec_inc_state_t *inc,
                        int                  S_orig,
                        const int           *prune_map,    /* nullable; NULL = identity */
                        const double        *in,
                        double              *gv,           /* length inc->n_indiv */
                        tree_pre_visit_fn    pre_visit,    /* nullable */
                        tree_visit_fn        visit,
                        void                *ctx,
                        int                  use_span);    /* 1 = weight G_s by span_s
                                                              (span-weighted eGRM, uniform
                                                              arms); 0 = plain G_s (learned
                                                              arms — span is non-identifiable
                                                              with free ws2 and ill-conditions
                                                              the CG). st.span = 1.0 when 0. */
