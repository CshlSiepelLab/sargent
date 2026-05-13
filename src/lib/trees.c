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
    int             vg_pos;        /* sites consumed so far (next call returns site vg_pos) */
    int            *poly_site_idx; /* [n_variants] poly index -> site index */
    int             cur_poly;
} trees_impl_t;

/* Count non-ancestral alleles for one diploid individual.
   Returns NaN if either haplotype has missing data. */
static double decode_indiv(int32_t h1, int32_t h2)
{
    if (h1 == TSK_MISSING_DATA || h2 == TSK_MISSING_DATA) return NAN;
    return (double)((h1 > 0) + (h2 > 0));
}

static void decode_ts_variant(const tsk_variant_t *var, int n_indiv, double *out)
{
    const int32_t *g = var->genotypes;
    for (int i = 0; i < n_indiv; i++)
        out[i] = decode_indiv(g[2 * i], g[2 * i + 1]);
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

/* Consume and discard sites until vg_pos == target.
   After return, the next tsk_vargen_next call yields site `target`. */
static int skip_to_site(trees_impl_t *t, int target)
{
    tsk_variant_t *var;
    while (t->vg_pos < target) {
        int ret = tsk_vargen_next(&t->vg, &var);
        if (ret < 0) { fprintf(stderr, "trees: %s\n", tsk_strerror(ret)); return -1; }
        if (ret == 0) return 0;  /* unexpectedly hit EOF */
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

        if (var->num_alleles != 2) continue;  /* skip multi-allelic */

        decode_ts_variant(var, n, tmp);

        double sum  = 0.0;
        int   n_obs = 0;
        for (int i = 0; i < n; i++) {
            if (!isnan(tmp[i])) { sum += tmp[i]; n_obs++; }
        }
        if (n_obs == 0) continue;

        double mean  = sum / n_obs;
        double p_hat = mean / 2.0;
        double var_  = 2.0 * p_hat * (1.0 - p_hat);
        if (var_ == 0.0) continue;  /* monomorphic */

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

    return init_vargen(t);  /* rewind for streaming */

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

    decode_ts_variant(var, r->n_samples, out);
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
    /* fast-forward to just before target; O(var_idx * n_samples) */
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
    free(t);
    free(r->stats);
    free(r);
    /* ts is NOT freed — caller owns it */
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

    trees_impl_t *t = calloc(1, sizeof(*t));
    if (!t) return NULL;
    t->ts            = ts;
    t->n_haplotypes  = n_hap;
    t->n_sites_total = tsk_treeseq_get_num_sites(ts);

    geno_reader_t *r = calloc(1, sizeof(*r));
    if (!r) { free(t); return NULL; }
    r->ops        = &trees_ops;
    r->n_samples  = (int)(n_hap / 2);
    r->n_variants = (int)t->n_sites_total;  /* updated to poly count by precompute_stats */
    r->impl       = t;
    return r;
}
