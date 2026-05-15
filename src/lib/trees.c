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

#include <tskit.h>
#include "geno.h"

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

    tsk_size_t nnd = tsk_treeseq_get_num_nodes(ts);
    double   *work     = malloc(2 * nnd * sizeof(double));
    tsk_id_t *order_buf = malloc(nnd * sizeof(tsk_id_t));
    double   *eps      = malloc((size_t)n_indiv * sizeof(double));
    double   *gv       = malloc((size_t)n_indiv * sizeof(double));
    if (!work || !order_buf || !eps || !gv) {
        free(work); free(order_buf); free(eps); free(gv);
        return 1.0;
    }

    /* Rademacher probe: each entry ±1 */
    unsigned int rng = (rng_seed != 0) ? rng_seed : 0xDEADBEEFu;
    for (int i = 0; i < n_indiv; i++) {
        rng ^= rng << 13; rng ^= rng >> 17; rng ^= rng << 5;
        eps[i] = (rng & 1u) ? 1.0 : -1.0;
    }

    double sum = 0.0;
    int    count = 0;
    tsk_tree_t tree;
    int ret = tsk_tree_init(&tree, ts, 0);
    if (ret < 0) goto cleanup;

    for (ret = tsk_tree_first(&tree); ret == TSK_TREE_OK;
         ret = tsk_tree_next(&tree)) {
        tree_cov_matvec(gv, eps, n_indiv, &tree, work, order_buf, hap_to_indiv);
        double dot = 0.0;
        for (int i = 0; i < n_indiv; i++) dot += eps[i] * gv[i];
        sum += dot / n_indiv;
        count++;
    }
    tsk_tree_free(&tree);

cleanup:
    free(work); free(order_buf); free(eps); free(gv);
    return (count > 0 && sum > 0.0) ? sum / count : 1.0;
}
