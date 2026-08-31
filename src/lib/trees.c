/*
 * SARGENT: Sparse ARG-based ENcoding of Traits
 *
 * Copyright (c) 2026, Cold Spring Harbor Laboratory
 * All rights reserved.
 *
 * This file is part of SARGENT and is distributed under the BSD 3-Clause License.
 * See the LICENSE file in the project root for details.
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <time.h>

#include <tskit.h>
#include "geno.h"
#include "sargent_omp.h"

#define TSK_CHECK(ret, label) do { \
    if ((ret) < 0) { \
        fprintf(stderr, "tskit error: %s\n", tsk_strerror(ret)); \
        goto label; \
    } \
} while (0)

typedef struct {
    tsk_treeseq_t  *ts;            /* borrowed — caller owns and frees */
    tsk_size_t      n_haplotypes;
    tsk_size_t      n_sites_total;
    tsk_vargen_t    vg;
    int             vg_ready;
    int             vg_pos;        /* sites consumed so far */
    int            *poly_site_idx; /* [n_variants] poly index -> site index */
    int             cur_poly;
    int            *hap_to_indiv;  /* [n_haplotypes] sample_rank -> individual */
    int            *indiv_haps;    /* [2*n_indiv] pairs of sample ranks per individual */
} trees_impl_t;

/* -------------------------------------------------------------------------
 * tsk_build_hap_to_indiv: build sample-rank → individual-index mapping.
 *
 * Uses the tree sequence's individual table when present (correct for any
 * sample ordering).  Falls back to rank/2 with a warning when no individual
 * table exists (valid only if msprime-style consecutive haplotype pairs are
 * guaranteed).
 * ---------------------------------------------------------------------- */
int *tsk_build_hap_to_indiv(const tsk_treeseq_t *ts)
{
    tsk_size_t n_hap   = tsk_treeseq_get_num_samples(ts);
    tsk_size_t n_indiv = tsk_treeseq_get_num_individuals(ts);

    int *h2i = malloc(n_hap * sizeof(int));
    if (!h2i) return NULL;

    if (n_indiv == 0) {
        /* No individual table — fall back to rank/2 */
        if (n_hap > 2)
            fprintf(stderr,
                    "trees: no individual table in tree sequence; "
                    "assuming consecutive rank/2 diploid pairing\n");
        for (tsk_size_t h = 0; h < n_hap; h++)
            h2i[h] = (int)(h / 2);
        return h2i;
    }

    if (n_indiv * 2 != n_hap) {
        fprintf(stderr,
                "trees: %lld individuals but %lld haplotypes; "
                "expected exactly 2 sample nodes per individual\n",
                (long long)n_indiv, (long long)n_hap);
        free(h2i);
        return NULL;
    }

    const tsk_id_t *sample_index_map = tsk_treeseq_get_sample_index_map(ts);

    for (tsk_size_t h = 0; h < n_hap; h++) h2i[h] = -1;

    for (tsk_size_t i = 0; i < n_indiv; i++) {
        tsk_size_t n_nodes = ts->individual_nodes_length[i];

        if (n_nodes != 2) {
            fprintf(stderr,
                    "trees: individual %lld has %lld nodes; expected 2 (diploid)\n",
                    (long long)i, (long long)n_nodes);
            free(h2i);
            return NULL;
        }

        for (tsk_size_t j = 0; j < n_nodes; j++) {
            tsk_id_t node = ts->individual_nodes[i][j];
            tsk_id_t rank = sample_index_map[node];
            if (rank < 0) {
                fprintf(stderr,
                        "trees: node %lld for individual %lld is not a sample\n",
                        (long long)node, (long long)i);
                free(h2i);
                return NULL;
            }
            h2i[(tsk_size_t)rank] = (int)i;
        }
    }

    return h2i;
}

/* Build indiv_haps[2*i], indiv_haps[2*i+1] = the two sample ranks for
 * individual i.  Derived by inverting hap_to_indiv. */
static int *build_indiv_haps(int n_hap, int n_indiv, const int *hap_to_indiv)
{
    int *ih = malloc((size_t)(2 * n_indiv) * sizeof(int));
    if (!ih) return NULL;
    for (int i = 0; i < 2 * n_indiv; i++) ih[i] = -1;

    int *slot = calloc((size_t)n_indiv, sizeof(int));
    if (!slot) { free(ih); return NULL; }

    for (int h = 0; h < n_hap; h++) {
        int ind = hap_to_indiv[h];
        if (ind < 0 || ind >= n_indiv) continue;
        int s = slot[ind]++;
        if (s < 2) ih[2 * ind + s] = h;
    }

    free(slot);
    return ih;
}

/* -------------------------------------------------------------------------
 * Genotype reader backend — tree sequence
 * ---------------------------------------------------------------------- */

/* Count non-ancestral alleles for one diploid individual.
   Returns NaN if either haplotype has missing data. */
static double decode_indiv(int32_t h1, int32_t h2)
{
    if (h1 == TSK_MISSING_DATA || h2 == TSK_MISSING_DATA) return NAN;
    return (double)((h1 > 0) + (h2 > 0));
}

/* Decode one variant into per-individual dosage using the explicit haplotype
 * mapping rather than the rank/2 assumption. */
static void decode_ts_variant(const tsk_variant_t *var, int n_indiv,
                               const int *indiv_haps, double *out)
{
    const int32_t *g = var->genotypes;
    for (int i = 0; i < n_indiv; i++)
        out[i] = decode_indiv(g[indiv_haps[2 * i]], g[indiv_haps[2 * i + 1]]);
}

static void apply_flags(double *out, int n, const geno_var_stats_t *s, int flags)
{
    int impute = flags & (GENO_IMPUTE | GENO_STANDARDIZE);
    int scale  = flags & GENO_STANDARDIZE;
    for (int i = 0; i < n; i++) {
        double g = out[i];
        if (isnan(g)) g = impute ? s->mean : NAN;
        if (scale)    g = (g - s->mean) * s->scale;
        out[i] = g;
    }
}

static int init_vargen(trees_impl_t *t)
{
    if (t->vg_ready) {
        tsk_vargen_free(&t->vg);
        t->vg_ready = 0;
    }
    const tsk_id_t *samples  = tsk_treeseq_get_samples(t->ts);
    tsk_size_t      n_samples = tsk_treeseq_get_num_samples(t->ts);
    int ret = tsk_vargen_init(&t->vg, t->ts, samples, n_samples, NULL, 0);
    if (ret < 0) {
        fprintf(stderr, "trees: tsk_vargen_init: %s\n", tsk_strerror(ret));
        return -1;
    }
    t->vg_ready = 1;
    t->vg_pos   = 0;
    return 0;
}

static int skip_to_site(trees_impl_t *t, int target)
{
    tsk_variant_t *var;
    while (t->vg_pos < target) {
        int ret = tsk_vargen_next(&t->vg, &var);
        if (ret < 0) { fprintf(stderr, "trees: %s\n", tsk_strerror(ret)); return -1; }
        if (ret == 0) return 0;
        t->vg_pos++;
    }
    return 1;
}

static int trees_precompute_stats(geno_reader_t *r)
{
    trees_impl_t *t     = r->impl;
    int           n     = r->n_samples;
    int           nsites = (int)t->n_sites_total;

    geno_var_stats_t *stats    = malloc((size_t)nsites * sizeof(*stats));
    int              *poly_idx = malloc((size_t)nsites * sizeof(int));
    double           *tmp      = malloc((size_t)n      * sizeof(double));
    if (!stats || !poly_idx || !tmp) goto oom;

    if (init_vargen(t) < 0) goto oom;

    tsk_variant_t *var;
    int n_poly = 0, ret;
    while ((ret = tsk_vargen_next(&t->vg, &var)) == 1) {
        int site_idx = t->vg_pos++;

        if (var->num_alleles != 2) continue;

        decode_ts_variant(var, n, t->indiv_haps, tmp);

        double sum  = 0.0;
        int   n_obs = 0;
        for (int i = 0; i < n; i++) {
            if (!isnan(tmp[i])) { sum += tmp[i]; n_obs++; }
        }
        if (n_obs == 0) continue;

        double mean  = sum / n_obs;
        double p_hat = mean / 2.0;
        double var_  = 2.0 * p_hat * (1.0 - p_hat);
        if (var_ == 0.0) continue;

        stats[n_poly].mean   = mean;
        stats[n_poly].scale  = 1.0 / sqrt(var_);
        stats[n_poly].n_miss = n - n_obs;
        poly_idx[n_poly]     = site_idx;
        n_poly++;
    }
    if (ret < 0) { fprintf(stderr, "trees: %s\n", tsk_strerror(ret)); goto oom; }

    free(tmp);
    free(r->stats);
    free(t->poly_site_idx);

    r->stats         = realloc(stats,    (size_t)n_poly * sizeof(*stats));
    t->poly_site_idx = realloc(poly_idx, (size_t)n_poly * sizeof(int));
    r->n_variants    = n_poly;
    t->cur_poly      = 0;

    return init_vargen(t);

oom:
    free(stats); free(poly_idx); free(tmp);
    return -1;
}

static int trees_read_variant(geno_reader_t *r, double *out, int flags)
{
    trees_impl_t *t = r->impl;
    if (!r->stats) {
        fprintf(stderr, "trees: call precompute_stats before reading\n");
        return -1;
    }
    if (t->cur_poly >= r->n_variants) return 0;

    int target = t->poly_site_idx[t->cur_poly];

    int ret = skip_to_site(t, target);
    if (ret <= 0) return ret;

    tsk_variant_t *var;
    ret = tsk_vargen_next(&t->vg, &var);
    if (ret < 0) { fprintf(stderr, "trees: %s\n", tsk_strerror(ret)); return -1; }
    if (ret == 0) return 0;
    t->vg_pos++;

    decode_ts_variant(var, r->n_samples, t->indiv_haps, out);
    if (flags) apply_flags(out, r->n_samples, &r->stats[t->cur_poly], flags);
    t->cur_poly++;
    return 1;
}

static int trees_read_chunk(geno_reader_t *r, int k, double *out, int flags)
{
    int n = r->n_samples, nread = 0;
    while (nread < k) {
        int ret = trees_read_variant(r, out + (size_t)nread * n, flags);
        if (ret <= 0) break;
        nread++;
    }
    return nread;
}

static int trees_seek(geno_reader_t *r, int var_idx)
{
    trees_impl_t *t = r->impl;
    if (var_idx < 0 || var_idx >= r->n_variants) return -1;
    if (init_vargen(t) < 0) return -1;
    t->cur_poly = var_idx;
    return skip_to_site(t, t->poly_site_idx[var_idx]) > 0 ? 0 : -1;
}

static void trees_rewind(geno_reader_t *r)
{
    trees_impl_t *t = r->impl;
    init_vargen(t);
    t->cur_poly = 0;
}

static void trees_close(geno_reader_t *r)
{
    trees_impl_t *t = r->impl;
    if (t->vg_ready) tsk_vargen_free(&t->vg);
    free(t->poly_site_idx);
    free(t->hap_to_indiv);
    free(t->indiv_haps);
    free(t);
    free(r->stats);
    free(r);
}

static const geno_ops_t trees_ops = {
    trees_precompute_stats,
    trees_read_variant,
    trees_read_chunk,
    trees_seek,
    trees_rewind,
    trees_close,
};

geno_reader_t *geno_open_trees(tsk_treeseq_t *ts)
{
    tsk_size_t n_hap = tsk_treeseq_get_num_samples(ts);
    if (n_hap % 2 != 0) {
        fprintf(stderr, "trees: odd number of haplotypes (%lld); diploid assumed\n",
                (long long)n_hap);
        return NULL;
    }
    int n_indiv = (int)(n_hap / 2);

    int *h2i = tsk_build_hap_to_indiv(ts);
    if (!h2i) return NULL;

    int *ih = build_indiv_haps((int)n_hap, n_indiv, h2i);
    if (!ih) { free(h2i); return NULL; }

    trees_impl_t *t = calloc(1, sizeof(*t));
    if (!t) { free(h2i); free(ih); return NULL; }
    t->ts            = ts;
    t->n_haplotypes  = n_hap;
    t->n_sites_total = tsk_treeseq_get_num_sites(ts);
    t->hap_to_indiv  = h2i;
    t->indiv_haps    = ih;

    geno_reader_t *r = calloc(1, sizeof(*r));
    if (!r) { free(h2i); free(ih); free(t); return NULL; }
    r->ops        = &trees_ops;
    r->n_samples  = n_indiv;
    r->n_variants = (int)t->n_sites_total;
    r->impl       = t;
    return r;
}

/* ---------------------------------------------------------------------------
 * tree_cov_matvec: compute out = G_s * in for a single local tree.
 *
 * G_s[i,j] = shared branch length between diploid individuals i and j
 *           = sum of branch lengths from root down to MRCA(i,j).
 *
 * Algorithm (O(n_nodes) two-pass tree traversal):
 *
 *   Post-order pass — subtree sums S[u]:
 *     leaf (sample haplotype h): S[u] = 0.5 * in[hap_to_indiv[rank(h)]]
 *       (factor 0.5 splits each individual's contribution across its two
 *        haplotype nodes; the pair sums to in[indiv] in the subtree)
 *     internal node u:  S[u] = sum of S[c] over children c
 *
 *   Pre-order pass — accumulated weighted path sums acc[u]:
 *     root r:  acc[r] = 0
 *     non-root u:  acc[u] = acc[parent[u]] + branch_length(u) * S[u]
 *
 *   Result: out[i] = acc[hap_a(i)] + acc[hap_b(i)]
 *
 * hap_to_indiv maps sample rank → individual index; build once with
 * tsk_build_hap_to_indiv and reuse across all tree iterations.
 * ---------------------------------------------------------------------------*/
void tree_cov_matvec(double *out, const double *in, int n_indiv,
                     const tsk_tree_t *tree,
                     double *work, tsk_id_t *order_buf,
                     const int *hap_to_indiv)
{
    const tsk_treeseq_t *ts  = tree->tree_sequence;
    tsk_size_t           nnd = tsk_treeseq_get_num_nodes(ts);

    double *S   = work;
    double *acc = work + nnd;

    const double   *time             = ts->tables->nodes.time;
    const tsk_id_t *sample_index_map = tsk_treeseq_get_sample_index_map(ts);
    const tsk_id_t *samples          = tsk_treeseq_get_samples(ts);

    for (tsk_size_t i = 0; i < nnd; i++) { S[i] = 0.0; acc[i] = 0.0; }

    /* Post-order: S[u] = weighted sum of in-vector over subtree of u */
    tsk_size_t n_order = 0;
    tsk_tree_postorder(tree, order_buf, &n_order);

    for (tsk_size_t k = 0; k < n_order; k++) {
        tsk_id_t u  = order_buf[k];
        tsk_id_t si = sample_index_map[u];
        if (si >= 0) {
            int indiv = hap_to_indiv[(int)si];
            if (indiv >= 0 && indiv < n_indiv)
                S[u] = 0.5 * in[indiv];
        } else {
            for (tsk_id_t c = tree->left_child[u]; c != TSK_NULL;
                 c = tree->right_sib[c])
                S[u] += S[c];
        }
    }

    /* Pre-order: acc[u] = cumulative (branch_length × S) from root to u */
    n_order = 0;
    tsk_tree_preorder(tree, order_buf, &n_order);

    for (tsk_size_t k = 0; k < n_order; k++) {
        tsk_id_t u = order_buf[k];
        tsk_id_t p = tree->parent[u];
        if (p == TSK_NULL) {
            acc[u] = 0.0;
        } else {
            double bl = time[p] - time[u];
            acc[u] = acc[p] + bl * S[u];
        }
    }

    /* Read off: out[i] = acc[hap_a(i)] + acc[hap_b(i)] */
    for (int i = 0; i < n_indiv; i++) out[i] = 0.0;

    tsk_size_t n_hap = tsk_treeseq_get_num_samples(ts);
    for (tsk_size_t h = 0; h < n_hap; h++) {
        int indiv = hap_to_indiv[(int)h];
        if (indiv >= 0 && indiv < n_indiv)
            out[indiv] += acc[samples[h]];
    }
}

/* ---------------------------------------------------------------------------
 * tree_cov_matvec_pairwise: REFERENCE implementation of out = G_s * in.
 *
 * Walks every haploid sample pair (h_a, h_b), computes MRCA(h_a, h_b) by
 * mark-and-walk, and accumulates the 0.5 * (time(root) - time(mrca)) entry
 * into the matvec result.  O(n_haps^2 * depth) per tree.
 *
 * Independent of the production two-pass post-/pre-order algorithm; used as
 * a paranoia cross-check (test/lib/check_tree_matvec) before we trust the
 * incremental edge-diff implementation in chunks 2-3.
 *
 * `mark` is caller-allocated scratch of length tsk_treeseq_get_num_nodes(ts);
 * contents are arbitrary on entry, clobbered on exit.
 * --------------------------------------------------------------------------- */
void tree_cov_matvec_pairwise(double *out, const double *in, int n_indiv,
                              const tsk_tree_t *tree,
                              int *mark,
                              const int *hap_to_indiv)
{
    const tsk_treeseq_t *ts      = tree->tree_sequence;
    const tsk_size_t     n_nodes = tsk_treeseq_get_num_nodes(ts);
    const tsk_id_t      *parent  = tree->parent;
    const double        *time    = ts->tables->nodes.time;
    const tsk_id_t      *samples = tsk_treeseq_get_samples(ts);
    const tsk_size_t     n_haps  = tsk_treeseq_get_num_samples(ts);

    for (int i = 0; i < n_indiv; i++) out[i] = 0.0;
    for (tsk_size_t u = 0; u < n_nodes; u++) mark[u] = -1;

    int version = 0;
    for (tsk_size_t a = 0; a < n_haps; a++) {
        tsk_id_t ha = samples[a];
        int ia = hap_to_indiv[(int)a];
        if (ia < 0 || ia >= n_indiv) continue;

        version++;
        /* Walk ha to root, marking ancestors and tracking the root. */
        tsk_id_t root_a = ha;
        for (tsk_id_t u = ha; u != TSK_NULL; u = parent[u]) {
            mark[u]  = version;
            root_a   = u;
        }
        double time_root = time[root_a];

        for (tsk_size_t b = 0; b < n_haps; b++) {
            tsk_id_t hb = samples[b];
            int ib = hap_to_indiv[(int)b];
            if (ib < 0 || ib >= n_indiv) continue;

            /* Walk hb to first marked node (= MRCA(ha, hb) within this tree). */
            tsk_id_t mrca = TSK_NULL;
            for (tsk_id_t u = hb; u != TSK_NULL; u = parent[u]) {
                if (mark[u] == version) { mrca = u; break; }
            }
            if (mrca == TSK_NULL) continue;  /* different roots (forest) */

            double T = time_root - time[mrca];
            out[ia] += 0.5 * T * in[ib];
        }
    }
}

/* ---------------------------------------------------------------------------
 * Incremental tree-cov matvec: maintains S[u], our_parent[u], bl[u] across
 * tree boundaries.  Phase 5 chunk 3a: no per-tree pre-order pass; matvec
 * read-out is a per-sample walk-up from sample to root, accumulating
 * bl[u] * S[u] along the way.
 *
 * On each tree boundary we read the edge diff from tree->tree_pos:
 *
 *   For each edge in OUT (parent p, child c) -- being removed:
 *     S[c] is unchanged (subtree below c is intact; only its attachment moves).
 *     Walk up our_parent[] from p to root, subtracting S[c] from each
 *     ancestor's S.  Then our_parent[c] := TSK_NULL, bl[c] := 0.
 *
 *   For each edge in IN (parent p, child c) -- being added:
 *     Set our_parent[c] := p, bl[c] := time[p] - time[c].
 *     Walk up our_parent[] from p to root, adding S[c] to each ancestor's S.
 *
 * Tskit's edge insertion/removal order guarantees consistency: at the time
 * an edge (p, c) is processed for removal our_parent[c] == p, and at the
 * time an edge (p, c) is processed for insertion the parent p is already
 * stitched into our_parent[].
 *
 * Matvec read-out:
 *   for each sample h:
 *       a = 0
 *       walk u = h, h's parent, ..., until our_parent[u] == TSK_NULL
 *           a += bl[u] * S[u]
 *       out[indiv(h)] += a
 *
 * O(n_haps * avg_depth) per tree; no n_nodes scan, no tsk_tree_preorder call,
 * no acc[] array.
 * --------------------------------------------------------------------------- */
/* Sort key for the one-time node ordering built in
 * tree_cov_matvec_inc_state_init: decreasing node time, id as a stable
 * tiebreak.  Ties can only occur between nodes with no edge between them
 * (tskit forbids time[parent] <= time[child]), so any tiebreak is safe. */
typedef struct { double t; tsk_id_t id; } node_time_key_t;

static int node_time_desc_cmp(const void *a, const void *b)
{
    const node_time_key_t *x = a, *y = b;
    if (x->t > y->t) return -1;
    if (x->t < y->t) return  1;
    return (x->id < y->id) ? -1 : (x->id > y->id);
}

static void build_tree_nodes(tree_cov_matvec_inc_state_t *st, int s);

int tree_cov_matvec_inc_state_init(tree_cov_matvec_inc_state_t *st,
                                    const tsk_treeseq_t *ts,
                                    int n_indiv,
                                    const int *hap_to_indiv)
{
    /* Zero EVERYTHING first.  The optional per-tree caches below are allocated
     * only when `compact` is set, but tree_cov_matvec_inc_state_free() frees
     * them unconditionally.  build_kinship_blocks memsets its state before
     * calling here, so it survived; gall_build declares its state on the stack
     * and does not, so those pointers held stack garbage and free() segfaulted
     * on the first cold .gall build.  An init function must fully initialise
     * its object rather than depend on the caller having zeroed it. */
    memset(st, 0, sizeof *st);

    st->ts            = ts;
    st->n_indiv       = n_indiv;
    st->hap_to_indiv  = hap_to_indiv;
    st->n_nodes       = tsk_treeseq_get_num_nodes(ts);

    st->S          = malloc((size_t)st->n_nodes * sizeof(double));
    st->our_parent = malloc((size_t)st->n_nodes * sizeof(tsk_id_t));
    st->bl         = malloc((size_t)st->n_nodes * sizeof(double));
    st->time_order = malloc((size_t)st->n_nodes * sizeof(tsk_id_t));
    st->cum        = malloc((size_t)st->n_nodes * sizeof(double));
    st->S0         = malloc((size_t)st->n_nodes * sizeof(double));
    if (!st->S || !st->our_parent || !st->bl || !st->time_order || !st->cum ||
        !st->S0) {
        tree_cov_matvec_inc_state_free(st);
        return -1;
    }

    /* Build the parent-before-child order once, for the whole tree sequence.
     * tskit enforces time[parent] > time[child] on every edge, so ordering
     * nodes by decreasing time is a valid topological order in every tree. */
    {
        const double *ntime = ts->tables->nodes.time;
        node_time_key_t *keys = malloc((size_t)st->n_nodes * sizeof(*keys));
        if (!keys) { tree_cov_matvec_inc_state_free(st); return -1; }
        for (tsk_size_t u = 0; u < st->n_nodes; u++) {
            keys[u].t  = ntime[u];
            keys[u].id = (tsk_id_t)u;
        }
        qsort(keys, (size_t)st->n_nodes, sizeof(*keys), node_time_desc_cmp);
        for (tsk_size_t k = 0; k < st->n_nodes; k++) st->time_order[k] = keys[k].id;
        free(keys);
    }

    /* Compact only when it pays.  A tree over n samples holds at most ~2n-1
     * nodes, so compare the node table against that.  Simulated ARG: 7,444 vs
     * ~6,000, ratio 1.2 -> sweep the global order, no cache, no extra memory.
     * Relate ARG: 368,878 vs ~6,000, ratio 61 -> compact, or the sweep does 60x
     * redundant work.  The threshold is deliberately loose; between 1 and 4 the
     * two paths cost about the same, so the choice cannot go badly wrong. */
    {
        tsk_size_t n_hap = tsk_treeseq_get_num_samples(ts);
        tsk_size_t per_tree = 2 * n_hap;
        st->compact = (per_tree > 0 && st->n_nodes > 4 * per_tree);
        st->cur_tree = 0;
        st->n_trees = tsk_treeseq_get_num_trees(ts);

        if (st->compact) {
            st->tree_nodes  = calloc((size_t)st->n_trees, sizeof(*st->tree_nodes));
            st->tree_nnodes = calloc((size_t)st->n_trees, sizeof(*st->tree_nnodes));
            st->has_child   = malloc((size_t)st->n_nodes * sizeof(*st->has_child));
            st->tree_parent = calloc((size_t)st->n_trees, sizeof(*st->tree_parent));
            st->tree_bl     = calloc((size_t)st->n_trees, sizeof(*st->tree_bl));
            st->tree_span   = calloc((size_t)st->n_trees, sizeof(*st->tree_span));
            st->struct_cached = 0;
            if (!st->tree_nodes || !st->tree_nnodes || !st->has_child ||
                !st->tree_parent || !st->tree_bl || !st->tree_span) {
                tree_cov_matvec_inc_state_free(st);
                return -1;
            }
        }
    }
    return 0;
}

void tree_cov_matvec_inc_state_free(tree_cov_matvec_inc_state_t *st)
{
    if (!st) return;
    free(st->S);          st->S          = NULL;
    free(st->our_parent); st->our_parent = NULL;
    free(st->bl);         st->bl         = NULL;
    free(st->time_order); st->time_order = NULL;
    free(st->cum);        st->cum        = NULL;
    free(st->S0);         st->S0         = NULL;
    if (st->tree_nodes) {
        for (tsk_size_t t = 0; t < st->n_trees; t++) free(st->tree_nodes[t]);
        free(st->tree_nodes); st->tree_nodes = NULL;
    }
    if (st->tree_parent) {
        for (tsk_size_t t = 0; t < st->n_trees; t++) free(st->tree_parent[t]);
        free(st->tree_parent); st->tree_parent = NULL;
    }
    if (st->tree_bl) {
        for (tsk_size_t t = 0; t < st->n_trees; t++) free(st->tree_bl[t]);
        free(st->tree_bl); st->tree_bl = NULL;
    }
    free(st->tree_span); st->tree_span = NULL;
    free(st->tree_nnodes); st->tree_nnodes = NULL;
    free(st->has_child);   st->has_child   = NULL;
}

void tree_cov_matvec_inc_state_reset(tree_cov_matvec_inc_state_t *st,
                                     const double *in)
{
    /* Tree index for the per-tree node cache.  The STATE counts its own
     * advances rather than trusting the caller to pass an index: the first
     * version required tree_walker_run_inc to set it, and check_tree_matvec --
     * which drives advance/readout directly -- left it pinned at 0, so every
     * tree reused tree 0's node list.  That is a silent wrong answer, and it
     * only surfaced because the validator happened to cover it. */
    st->cur_tree = -1;

    const tsk_id_t *sample_index_map = tsk_treeseq_get_sample_index_map(st->ts);

    for (tsk_size_t u = 0; u < st->n_nodes; u++) {
        st->our_parent[u] = TSK_NULL;
        st->bl[u]         = 0.0;
        tsk_id_t si = sample_index_map[u];
        if (si >= 0) {
            int indiv = st->hap_to_indiv[(int)si];
            if (indiv >= 0 && indiv < st->n_indiv) {
                st->S[u] = 0.5 * in[indiv];
            } else {
                st->S[u] = 0.0;
            }
        } else {
            st->S[u] = 0.0;
        }
        st->S0[u] = st->S[u];      /* seed kept for the replay path */
    }
}

/* Apply the edge diff at tree->tree_pos to S[], our_parent[], bl[].
 * Must be called for every tree boundary even if the readout is skipped --
 * S[] tracks the post-order subtree sum and goes out of sync otherwise. */
void tree_cov_matvec_inc_advance(tree_cov_matvec_inc_state_t *st,
                                  const tsk_tree_t *tree)
{
    /* Tree index for the per-tree node cache.  The state counts its own
     * advances rather than trusting a caller to pass an index: requiring that
     * left check_tree_matvec -- which drives advance/readout directly -- pinned
     * at tree 0, so every tree reused tree 0's node list.  A silent wrong
     * answer, caught only because the validator covers that path. */
    st->cur_tree++;

    const tsk_table_collection_t *tables      = st->ts->tables;
    const tsk_id_t               *edge_parent = tables->edges.parent;
    const tsk_id_t               *edge_child  = tables->edges.child;
    const double                 *time        = tables->nodes.time;
    const tsk_tree_position_t    *pos         = &tree->tree_pos;

    for (tsk_id_t j = pos->out.start; j != pos->out.stop; j++) {
        tsk_id_t e = pos->out.order[j];
        tsk_id_t p = edge_parent[e];
        tsk_id_t c = edge_child[e];
        double   sc = st->S[c];
        for (tsk_id_t u = p; u != TSK_NULL; u = st->our_parent[u])
            st->S[u] -= sc;
        st->our_parent[c] = TSK_NULL;
        st->bl[c]         = 0.0;
    }
    for (tsk_id_t j = pos->in.start; j != pos->in.stop; j++) {
        tsk_id_t e = pos->in.order[j];
        tsk_id_t p = edge_parent[e];
        tsk_id_t c = edge_child[e];
        st->our_parent[c] = p;
        st->bl[c]         = time[p] - time[c];
        double sc = st->S[c];
        for (tsk_id_t u = p; u != TSK_NULL; u = st->our_parent[u])
            st->S[u] += sc;
    }
}

/* Per-sample walk-up: out[indiv(h)] += sum_{v ancestor of h} bl[v] * S[v].
 *
 * Uses #pragma omp for (worksharing) rather than #pragma omp parallel for
 * (fork/join) so it participates in an EXISTING parallel region from
 * tree_walker_run_inc without spawning a nested team.  When called outside
 * any parallel region (e.g. from tree_cov_matvec_incremental), the orphaned
 * omp for is a no-op and the loop runs serially — correct but unparallelised.
 *
 * Each individual has exactly 2 haplotypes so the omp atomic below is at
 * most 2-way contended — cheaper than a per-thread scratch array. */
/* Root-path prefix sum, shared by both readout entry points.
 *
 * The obvious readout walks every sample independently from its leaf to the
 * root.  That costs sum-of-sample-depths per tree and makes sibling samples
 * re-traverse the same shared ancestral path over and over: on a 3000-sample,
 * 7444-node ARG it is 52,379 steps per tree against 5,998 edges actually
 * present -- the same path arithmetic done ~8.7x.
 *
 * Instead sweep the nodes once in parent-before-child order and carry the
 * running root-path sum down:
 *
 *     cum[u] = bl[u]*S[u] + cum[parent(u)],   cum[root] = 0
 *
 * so cum[sample] is exactly the quantity the leaf-to-root loop produced.
 * `time_order` is built once at init (see geno.h) and is valid for every tree,
 * so this costs a flat n_nodes pass with no per-tree setup.
 *
 * NOTE the sum is now right-associated (root-end first) where the old loop was
 * left-associated (leaf-end first).  The value is mathematically identical but
 * the rounding differs, so this is not a bit-identical change; it is checked
 * against a tolerance instead, and measures MORE accurate than the loop it
 * replaced (each path sum is formed once, not re-derived per descendant). */
/* Time-ordered list of the nodes present in the current tree, built once per
 * tree index and reused across the n_tr column passes (they all traverse the
 * same sequence, so the structure repeats even though S[] does not).
 *
 * "Present" = has a parent, or is a root with at least one child, OR IS A
 * SAMPLE.  All three matter:
 *   - roots, because the sweep is what zeroes cum[root] and a child reads it;
 *   - samples unconditionally, because an inferred ARG can leave a sample
 *     ISOLATED in a given tree (no parent, no children).  The global sweep gave
 *     such a node cum = 0 for free; a compacted list that skips it leaves
 *     cum[sample] holding a STALE value from the previous tree, and the gather
 *     loop reads exactly that.  This was not hypothetical -- it produced
 *     naive=-2198.47 against inc=1.36 on a Relate cell. */
static void build_tree_nodes(tree_cov_matvec_inc_state_t *st, int s)
{
    const tsk_id_t *parent  = st->our_parent;
    unsigned char  *hc      = st->has_child;
    const tsk_id_t *samples = tsk_treeseq_get_samples(st->ts);
    const tsk_size_t n_hap  = tsk_treeseq_get_num_samples(st->ts);

    memset(hc, 0, (size_t)st->n_nodes);
    for (tsk_size_t u = 0; u < st->n_nodes; u++)
        if (parent[u] != TSK_NULL) hc[parent[u]] = 1;
    for (tsk_size_t hi = 0; hi < n_hap; hi++) hc[samples[hi]] = 1;

    tsk_size_t cnt = 0;
    for (tsk_size_t u = 0; u < st->n_nodes; u++)
        if (parent[u] != TSK_NULL || hc[u]) cnt++;

    tsk_id_t *list = malloc((size_t)(cnt ? cnt : 1) * sizeof(*list));
    if (!list) return;                    /* caller falls back to the global sweep */
    tsk_size_t m = 0;
    for (tsk_size_t k = 0; k < st->n_nodes; k++) {
        tsk_id_t u = st->time_order[k];
        if (parent[u] != TSK_NULL || hc[u]) list[m++] = u;
    }
    st->tree_nodes[s]  = list;
    st->tree_nnodes[s] = m;

    /* Structure replay: keep parent and branch length for exactly these nodes,
     * so later columns can restore the tree without touching tskit. */
    if (st->tree_parent && st->tree_bl) {
        tsk_id_t *pp = malloc((size_t)(m ? m : 1) * sizeof(*pp));
        double   *bb = malloc((size_t)(m ? m : 1) * sizeof(*bb));
        if (pp && bb) {
            for (tsk_size_t k = 0; k < m; k++) { pp[k] = parent[list[k]]; bb[k] = st->bl[list[k]]; }
            st->tree_parent[s] = pp;
            st->tree_bl[s]     = bb;
        } else { free(pp); free(bb); }
    }
}

/* Root-path prefix sum: cum[u] = bl[u]*S[u] + cum[parent(u)], cum[root] = 0,
 * swept in an order that puts every parent before its children. */
static inline void readout_prefix(const tree_cov_matvec_inc_state_t *cst)
{
    /* The lazily-built per-tree list is a cache, so this mutates through a
     * const handle.  The readout's observable output is unchanged. */
    tree_cov_matvec_inc_state_t *st = (tree_cov_matvec_inc_state_t *)cst;

    const tsk_id_t *ord = st->time_order;
    tsk_size_t      n   = st->n_nodes;

    if (st->compact && st->cur_tree >= 0 &&
        (tsk_size_t)st->cur_tree < st->n_trees) {
        int t = st->cur_tree;
        if (!st->tree_nodes[t]) build_tree_nodes(st, t);
        if (st->tree_nodes[t]) { ord = st->tree_nodes[t]; n = st->tree_nnodes[t]; }
    }

    const tsk_id_t *parent = st->our_parent;
    const double   *bl     = st->bl;
    const double   *S      = st->S;
    double         *cum    = st->cum;

    for (tsk_size_t k = 0; k < n; k++) {
        tsk_id_t u = ord[k];
        tsk_id_t p = parent[u];
        cum[u] = (p == TSK_NULL) ? 0.0 : cum[p] + bl[u] * S[u];
    }
}

/* Serial readout.  Contains NO OpenMP directives, deliberately.
 *
 * build_kinship_blocks parallelises over COLUMNS, so this runs with each
 * thread inside its own column.  Orphaned `omp for` / `omp single` here would
 * bind to that enclosing team, and the threads reach them while sitting in
 * DIFFERENT columns -- worksharing regions must be encountered by every thread
 * of a team in the same order, so that is undefined behaviour.  It corrupted
 * the heap in practice: correct at 1 thread, "prediction stage failed" at 2.
 * Use readout_team() when, and only when, the calling team is walking one
 * tree sequence together. */
void tree_cov_matvec_inc_readout(double *out,
                                  const tree_cov_matvec_inc_state_t *st)
{
    const int n_indiv = st->n_indiv;
    for (int i = 0; i < n_indiv; i++) out[i] = 0.0;

    readout_prefix(st);

    const tsk_id_t  *samples = tsk_treeseq_get_samples(st->ts);
    const tsk_size_t n_hap   = tsk_treeseq_get_num_samples(st->ts);
    for (tsk_size_t hi = 0; hi < n_hap; hi++) {
        int indiv = st->hap_to_indiv[(int)hi];
        if (indiv < 0 || indiv >= n_indiv) continue;
        out[indiv] += st->cum[samples[hi]];
    }
}

#ifdef SARGENT_HAS_OPENMP
/* Worksharing readout, for the walker's own hoisted team only: every thread of
 * that team is on the SAME tree, so it reaches these constructs in step.
 *
 * Zeroing must complete before any thread accumulates into out[], or a
 * straggler still zeroing out[k] clobbers a faster thread's `out[k] +=` --
 * that corrupted the readout and made G_all indefinite (issue #26).  The
 * `omp single` exit barrier supplies that ordering, and folding the zeroing
 * inside it keeps the barrier count at two per readout. */
static void readout_team(double *out, const tree_cov_matvec_inc_state_t *st)
{
    const int n_indiv = st->n_indiv;
    const tsk_id_t  *samples = tsk_treeseq_get_samples(st->ts);
    const tsk_size_t n_hap   = tsk_treeseq_get_num_samples(st->ts);

    #pragma omp single
    {
        for (int i = 0; i < n_indiv; i++) out[i] = 0.0;
        readout_prefix(st);
    }   /* implicit barrier: out[] zeroed and cum[] complete for all threads */

    #pragma omp for schedule(static)
    for (tsk_size_t hi = 0; hi < n_hap; hi++) {
        int indiv = st->hap_to_indiv[(int)hi];
        if (indiv < 0 || indiv >= n_indiv) continue;
        #pragma omp atomic
        out[indiv] += st->cum[samples[hi]];
    }
}
#endif

void tree_cov_matvec_incremental(double *out,
                                 tree_cov_matvec_inc_state_t *st,
                                 const tsk_tree_t *tree)
{
    tree_cov_matvec_inc_advance(st, tree);
    tree_cov_matvec_inc_readout(out, st);
}

/* Mean genomic span per tree = sequence_length / num_trees.  Used to normalise
 * the per-tree span weight so Sum_s (span_s / mean_span) == S_orig: the
 * span-weighted kernel Sum_s span_s G_s then keeps the same overall scale as the
 * old unweighted Sum_s G_s (only the RELATIVE per-tree weighting changes). */
static double tree_seq_mean_span(const tsk_treeseq_t *ts)
{
    double     L  = tsk_treeseq_get_sequence_length(ts);
    tsk_size_t nt = tsk_treeseq_get_num_trees(ts);
    return (nt > 0 && L > 0.0) ? L / (double)nt : 1.0;
}

/* ---------------------------------------------------------------------------
 * tree_walker_run_inc: incremental tree-walker driving the advance/readout
 * split.  pre_visit may skip the readout; advance always runs to keep S[]
 * coherent across the tree-sequence walk.
 * --------------------------------------------------------------------------- */
/* tree_walker_run_inc: incremental tree-walker with hoisted OpenMP parallel
 * region.
 *
 * Phase perf: the previous design called tree_cov_matvec_inc_readout (which
 * had #pragma omp parallel for inside) once per tree — spawning and joining
 * a thread team ~S times per matvec, ~700S times per EM iteration.  At
 * S=3000 that is >2M fork/join events per run; on Linux pthreads each costs
 * ~5 µs → ~10 s overhead alone, masking the actual computation.
 *
 * Fix: hoist #pragma omp parallel outside the tree loop so one team is
 * created per tree_walker_run_inc call.  Inside:
 *   - omp single (with implicit barrier) handles the sequential steps
 *     (tsk_tree_next, inc_advance, pre_visit, visit callback).
 *   - omp for (worksharing within the existing team) handles readout.
 *
 * The serial fallback (no OpenMP, or n_hap ≤ 1000) is unchanged. */
int tree_walker_run_inc(tree_cov_matvec_inc_state_t *inc,
                        int                  S_orig,
                        const int           *prune_map,
                        const double        *in,
                        double              *gv,
                        tree_pre_visit_fn    pre_visit,
                        tree_visit_fn        visit,
                        void                *ctx,
                        int                  use_span)
{
    tsk_tree_t tree;
    int ret = tsk_tree_init(&tree, inc->ts, 0);
    if (ret < 0) {
        fprintf(stderr, "tree_walker_run_inc: tsk_tree_init: %s\n", tsk_strerror(ret));
        return -1;
    }

    tree_cov_matvec_inc_state_reset(inc, in);

    tree_walk_state_t st = {
        .n_indiv = inc->n_indiv,
        .span    = 1.0,
        .gv      = gv,
    };

    /* Span weighting: aggregation visitors weight each tree's G_s by st.span
     * (= span_s / mean_span), making the kernel the span-weighted eGRM. */
    double mean_span = tree_seq_mean_span(inc->ts);

    int err = 0;

    /* ---- structure replay -------------------------------------------------
     * The tree sequence does not change between columns, but tsk_tree_next
     * re-derives it every time -- 77% of the kinship build on an inferred ARG.
     * Once the first traversal has cached every tree, later columns restore the
     * structure from that cache and never call tskit again.  S has to be
     * rebuilt from the leaves rather than by edge diffs, which is O(nodes in
     * tree) and cheaper than the diffs were. */
    if (inc->compact && inc->struct_cached && !prune_map &&
        S_orig == (int)inc->n_trees) {
        const double mspan = mean_span;
        for (int t = 0; t < S_orig; t++) {
            const tsk_id_t *ord = inc->tree_nodes[t];
            const tsk_id_t *par = inc->tree_parent[t];
            const double   *bl  = inc->tree_bl[t];
            tsk_size_t      n   = inc->tree_nnodes[t];
            if (!ord || !par || !bl) { err = 1; break; }

            for (tsk_size_t k = 0; k < n; k++) {
                tsk_id_t u = ord[k];
                inc->our_parent[u] = par[k];
                inc->bl[u]         = bl[k];
                inc->S[u]          = inc->S0[u];
            }
            /* bottom-up: ord is decreasing in time, so backwards visits every
             * child before its parent */
            for (tsk_size_t k = n; k-- > 0; ) {
                tsk_id_t u = ord[k], pu = par[k];
                if (pu != TSK_NULL) inc->S[pu] += inc->S[u];
            }
            inc->cur_tree = t;
            tree_cov_matvec_inc_readout(gv, inc);
            st.s_orig = t;
            st.sp     = t;
            st.span   = use_span ? inc->tree_span[t] / mspan : 1.0;
            visit(&st, ctx);
        }
        tsk_tree_free(&tree);
        return err ? -1 : 0;
    }

#ifdef SARGENT_HAS_OPENMP
    tsk_size_t n_hap = tsk_treeseq_get_num_samples(inc->ts);

    /* `omp_in_parallel()`: when build_kinship_blocks parallelises over COLUMNS,
     * each thread already owns a whole column and walks the sequence alone.
     * Forking a second team here would nest, so take the serial path instead --
     * the column loop is the better place to spend the threads, because the
     * readout's prefix sweep is a serial dependence (a node needs its parent)
     * and does not workshare. */
    if ((int)n_hap > 1000 && !SARGENT_OMP_IN_PARALLEL()) {
        /* Hoisted parallel region: one fork for the entire tree sequence.
         * Shared control variables written only inside omp single blocks;
         * all threads see consistent values after the implicit barrier. */
        int s_shared   = 0;
        int ret_shared = TSK_TREE_OK;
        int do_visit   = 0;    /* -1 = exit signal; 0 = skip; 1 = readout+visit */

        #pragma omp parallel shared(s_shared, ret_shared, do_visit, err, tree, st)
        {
            /* Position to first tree (single thread). */
            #pragma omp single
            { ret_shared = tsk_tree_first(&tree); }

            for (;;) {
                /* Sequential: advance state and decide whether to readout. */
                #pragma omp single
                {
                    if (ret_shared != TSK_TREE_OK || s_shared >= S_orig) {
                        do_visit = -1;   /* signal loop exit to all threads */
                    } else {
                        int sp = prune_map ? prune_map[s_shared] : s_shared;
                        tree_cov_matvec_inc_advance(inc, &tree);
                        st.s_orig = s_shared;
                        st.sp     = sp;
                        st.span   = use_span ? (tree.interval.right - tree.interval.left) / mean_span : 1.0;
                        do_visit  = (!pre_visit || pre_visit(s_shared, sp, ctx))
                                    ? 1 : 0;
                        s_shared++;
                        /* Advance tree for the next iteration. */
                        ret_shared = tsk_tree_next(&tree);
                    }
                }
                /* All threads join here before checking do_visit. */

                if (do_visit < 0) break;   /* all threads exit */

                if (do_visit > 0) {
                    /* Worksharing readout: every thread of this team is on
                     * the same tree, so readout_team's `omp for` is safe here. */
                    readout_team(gv, inc);
                    /* Sequential callback (may modify ctx). */
                    #pragma omp single
                    { visit(&st, ctx); }
                }
            }

            /* Check for real iteration error (not just normal end-of-sequence). */
            #pragma omp single
            {
                if (ret_shared < 0 && ret_shared != 0   /* end-of-sequence: tsk_tree_next returns 0 (not TSK_TREE_OK=1) */) {
                    fprintf(stderr, "tree_walker_run_inc: tskit error: %s\n",
                            tsk_strerror(ret_shared));
                    err = 1;
                }
            }
        }
        /* end parallel region */

    } else {
        /* Serial path for small n_hap (fork overhead > benefit). */
        int s = 0;
        for (ret = tsk_tree_first(&tree); ret == TSK_TREE_OK && s < S_orig;
             ret = tsk_tree_next(&tree), s++) {
            int sp = prune_map ? prune_map[s] : s;
            tree_cov_matvec_inc_advance(inc, &tree);
            if (pre_visit && !pre_visit(s, sp, ctx)) continue;
            tree_cov_matvec_inc_readout(gv, inc);
            st.s_orig = s; st.sp = sp;
            double raw_span = tree.interval.right - tree.interval.left;
            if (inc->compact && inc->tree_span && s < (int)inc->n_trees)
                inc->tree_span[s] = raw_span;
            st.span = use_span ? raw_span / mean_span : 1.0;
            visit(&st, ctx);
        }
        if (ret < 0) {
            fprintf(stderr, "tree_walker_run_inc: tskit error: %s\n",
                    tsk_strerror(ret));
            err = 1;
        }
        /* Arm the replay path only when this traversal saw the WHOLE sequence
         * and every tree got cached.  A pruned or truncated walk leaves holes,
         * and replaying those would silently use another tree's structure. */
        if (!err && inc->compact && !prune_map && s == (int)inc->n_trees &&
            S_orig >= (int)inc->n_trees) {
            int all = 1;
            for (tsk_size_t t = 0; t < inc->n_trees && all; t++)
                if (!inc->tree_nodes[t] || !inc->tree_parent[t] || !inc->tree_bl[t])
                    all = 0;
            inc->struct_cached = all;
        }
    }

#else
    /* No OpenMP: original serial implementation. */
    {
        int s = 0;
        for (ret = tsk_tree_first(&tree); ret == TSK_TREE_OK && s < S_orig;
             ret = tsk_tree_next(&tree), s++) {
            int sp = prune_map ? prune_map[s] : s;
            tree_cov_matvec_inc_advance(inc, &tree);
            if (pre_visit && !pre_visit(s, sp, ctx)) continue;
            tree_cov_matvec_inc_readout(gv, inc);
            st.s_orig = s; st.sp = sp;
            double raw_span = tree.interval.right - tree.interval.left;
            if (inc->compact && inc->tree_span && s < (int)inc->n_trees)
                inc->tree_span[s] = raw_span;
            st.span = use_span ? raw_span / mean_span : 1.0;
            visit(&st, ctx);
        }
        if (ret < 0) {
            fprintf(stderr, "tree_walker_run_inc: tskit error: %s\n",
                    tsk_strerror(ret));
            err = 1;
        }
        /* Arm the replay path only when this traversal saw the WHOLE sequence
         * and every tree got cached.  A pruned or truncated walk leaves holes,
         * and replaying those would silently use another tree's structure. */
        if (!err && inc->compact && !prune_map && s == (int)inc->n_trees &&
            S_orig >= (int)inc->n_trees) {
            int all = 1;
            for (tsk_size_t t = 0; t < inc->n_trees && all; t++)
                if (!inc->tree_nodes[t] || !inc->tree_parent[t] || !inc->tree_bl[t])
                    all = 0;
            inc->struct_cached = all;
        }
    }
#endif

    tsk_tree_free(&tree);
    return err ? -1 : 0;
}

/* ---------------------------------------------------------------------------
 * tree_walker_run: iterate a tree sequence, invoking a per-tree visitor.
 *
 * Phase 5 chunk 1: this implementation is a thin wrapper around the existing
 * tsk_tree_first/_next loop + tree_cov_matvec.  No behaviour change; the
 * value-add is centralising tree-iteration across the three callers so the
 * incremental tree-cov computation (chunks 2-3) can replace the inner work
 * without touching call sites.
 * --------------------------------------------------------------------------- */
int tree_walker_run(const tsk_treeseq_t *ts,
                    const int           *hap_to_indiv,
                    int                  n_indiv,
                    int                  S_orig,
                    const int           *prune_map,
                    const double        *in,
                    double              *gv,
                    double              *work,
                    tsk_id_t            *order_buf,
                    tree_pre_visit_fn    pre_visit,
                    tree_visit_fn        visit,
                    void                *ctx)
{
    tsk_tree_t tree;
    int ret = tsk_tree_init(&tree, ts, 0);
    if (ret < 0) {
        fprintf(stderr, "tree_walker_run: tsk_tree_init: %s\n", tsk_strerror(ret));
        return -1;
    }

    tree_walk_state_t st = {
        .n_indiv = n_indiv,
        .s_orig  = 0,
        .sp      = 0,
        .span    = 1.0,
        .gv      = gv,
    };

    double mean_span = tree_seq_mean_span(ts);

    int s = 0;
    for (ret = tsk_tree_first(&tree); ret == TSK_TREE_OK && s < S_orig;
         ret = tsk_tree_next(&tree), s++) {
        int sp = prune_map ? prune_map[s] : s;
        if (pre_visit && !pre_visit(s, sp, ctx))
            continue;
        tree_cov_matvec(gv, in, n_indiv, &tree, work, order_buf, hap_to_indiv);
        st.s_orig = s;
        st.sp     = sp;
        st.span   = (tree.interval.right - tree.interval.left) / mean_span;
        visit(&st, ctx);
    }
    if (ret < 0) {
        fprintf(stderr, "tree_walker_run: tsk_tree iteration: %s\n", tsk_strerror(ret));
        tsk_tree_free(&tree);
        return -1;
    }
    tsk_tree_free(&tree);
    return 0;
}

/* Context + visitor for estimate_g_scale below. */
typedef struct {
    const double *eps;
    int           n_indiv;
    double        sum;
    int           count;
} g_scale_ctx_t;

static void g_scale_visit(const tree_walk_state_t *st, void *vctx)
{
    g_scale_ctx_t *c = vctx;
    double dot = 0.0;
    for (int i = 0; i < c->n_indiv; i++) dot += c->eps[i] * st->gv[i];
    c->sum  += st->span * dot / c->n_indiv;
    c->count++;
}

/* Estimate E[tr(G_s)/n] using a single Rademacher probe over all trees.
 *
 * Performs one complete tree-sequence walk; for each tree s computes
 *   eps^T (G_s eps) / n
 * which is an unbiased estimator of tr(G_s)/n.  Averaging over all S trees
 * gives an estimate of the mean per-tree kinship scale.  This scale factor
 * is used to calibrate the tau2_mix prior grid so that it is expressed in
 * phenotype-variance units regardless of branch-length scale.
 *
 * Returns 1.0 on allocation failure or empty tree sequence. */
double estimate_g_scale(const tsk_treeseq_t *ts, const int *hap_to_indiv,
                        int n_indiv, unsigned int rng_seed)
{
    int S = (int)tsk_treeseq_get_num_trees(ts);
    if (S <= 0 || n_indiv <= 0) return 1.0;

    double *eps = malloc((size_t)n_indiv * sizeof(double));
    double *gv  = malloc((size_t)n_indiv * sizeof(double));
    tree_cov_matvec_inc_state_t inc;
    int inc_ok = (tree_cov_matvec_inc_state_init(&inc, ts, n_indiv, hap_to_indiv) == 0);
    if (!eps || !gv || !inc_ok) {
        free(eps); free(gv);
        if (inc_ok) tree_cov_matvec_inc_state_free(&inc);
        return 1.0;
    }

    /* Rademacher probe: each entry ±1 */
    unsigned int rng = (rng_seed != 0) ? rng_seed : 0xDEADBEEFu;
    for (int i = 0; i < n_indiv; i++) {
        rng ^= rng << 13; rng ^= rng >> 17; rng ^= rng << 5;
        eps[i] = (rng & 1u) ? 1.0 : -1.0;
    }

    g_scale_ctx_t gctx = { .eps = eps, .n_indiv = n_indiv, .sum = 0.0, .count = 0 };
    int ret = tree_walker_run_inc(&inc, S, /*prune_map=*/NULL,
                                  eps, gv, NULL, g_scale_visit, &gctx,
                                  /*use_span=*/0);  /* learned prior-grid scale: un-span */
    double sum   = gctx.sum;
    int    count = gctx.count;

    tree_cov_matvec_inc_state_free(&inc);
    free(eps); free(gv);
    (void)ret;
    return (count > 0 && sum > 0.0) ? sum / count : 1.0;
}
