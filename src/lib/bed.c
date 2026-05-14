/*
 * SARGENT: Sparse ARG-based ENcoding of Traits
 *
 * Copyright (c) 2026, Cold Spring Harbor Laboratory
 * All rights reserved.
 *
 * This file is part of SARGENT and is distributed under the BSD 3-Clause License.
 * See the LICENSE file in the project root for details.
 */

#define _FILE_OFFSET_BITS 64

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <errno.h>

#include "geno.h"

static const unsigned char BED_MAGIC[3] = {0x6c, 0x1b, 0x01};

/* BED 2-bit encoding: 00=hom-ref(0), 01=missing(NaN), 10=het(1), 11=hom-alt(2) */
static const double GENO_LUT[4] = {0.0, NAN, 1.0, 2.0};

typedef struct {
    FILE    *fp;
    int      n_samples_total;   /* from .fam */
    int      n_variants_total;  /* from .bim */
    size_t   bytes_per_var;     /* ceil(n_samples_total / 4) */
    uint8_t *raw_buf;           /* scratch buffer, bytes_per_var bytes */
    int     *poly_idx;          /* [n_variants] poly index -> file variant index */
    int      cur_poly;          /* streaming position in poly variant space */
} bed_impl_t;

static int count_lines(const char *path)
{
    FILE *f = fopen(path, "r");
    if (!f) return -1;
    int n = 0, last = '\n', c;
    while ((c = fgetc(f)) != EOF) {
        if (c == '\n') n++;
        last = c;
    }
    if (last != '\n') n++;  /* no trailing newline */
    fclose(f);
    return n;
}

static void decode_raw(const uint8_t *raw, int n, double *out)
{
    for (int i = 0; i < n; i++)
        out[i] = GENO_LUT[(raw[i / 4] >> (2 * (i % 4))) & 0x3];
}

static void apply_flags(double *out, int n, const geno_var_stats_t *s, int flags)
{
    int impute = flags & (GENO_IMPUTE | GENO_STANDARDIZE);
    int scale  = flags & GENO_STANDARDIZE;
    for (int i = 0; i < n; i++) {
        double g = out[i];
        if (isnan(g))  g = impute ? s->mean : NAN;
        if (scale)     g = (g - s->mean) * s->scale;
        out[i] = g;
    }
}

static int bed_precompute_stats(geno_reader_t *r)
{
    bed_impl_t *b = r->impl;
    int n = b->n_samples_total;
    int p = b->n_variants_total;

    geno_var_stats_t *stats   = malloc((size_t)p * sizeof(*stats));
    int              *poly_idx = malloc((size_t)p * sizeof(int));
    double           *tmp      = malloc((size_t)n * sizeof(double));
    if (!stats || !poly_idx || !tmp) goto oom;

    if (fseeko(b->fp, 3, SEEK_SET) != 0) goto oom;

    int n_poly = 0;
    for (int vi = 0; vi < p; vi++) {
        if (fread(b->raw_buf, 1, b->bytes_per_var, b->fp) != b->bytes_per_var)
            goto oom;

        decode_raw(b->raw_buf, n, tmp);

        double sum = 0.0;
        int n_obs = 0;
        for (int i = 0; i < n; i++) {
            if (!isnan(tmp[i])) { sum += tmp[i]; n_obs++; }
        }
        if (n_obs == 0) continue;

        double mean  = sum / n_obs;
        double p_hat = mean / 2.0;
        double var   = 2.0 * p_hat * (1.0 - p_hat);
        if (var == 0.0) continue;   /* monomorphic: skip */

        stats[n_poly].mean   = mean;
        stats[n_poly].scale  = 1.0 / sqrt(var);
        stats[n_poly].n_miss = n - n_obs;
        poly_idx[n_poly]     = vi;
        n_poly++;
    }

    free(tmp);
    free(r->stats);
    free(b->poly_idx);

    r->stats      = realloc(stats,    (size_t)n_poly * sizeof(*stats));
    b->poly_idx   = realloc(poly_idx, (size_t)n_poly * sizeof(int));
    r->n_variants = n_poly;
    b->cur_poly   = 0;

    if (n_poly > 0) {
        off_t off = 3 + (off_t)b->poly_idx[0] * (off_t)b->bytes_per_var;
        fseeko(b->fp, off, SEEK_SET);
    }
    return 0;

oom:
    free(stats); free(poly_idx); free(tmp);
    return -1;
}

static int bed_seek(geno_reader_t *r, int var_idx)
{
    bed_impl_t *b = r->impl;
    if (var_idx < 0 || var_idx >= r->n_variants) return -1;
    off_t off = 3 + (off_t)b->poly_idx[var_idx] * (off_t)b->bytes_per_var;
    if (fseeko(b->fp, off, SEEK_SET) != 0) return -1;
    b->cur_poly = var_idx;
    return 0;
}

static void bed_rewind(geno_reader_t *r)
{
    bed_seek(r, 0);
}

static int bed_read_variant(geno_reader_t *r, double *out, int flags)
{
    bed_impl_t *b = r->impl;
    if (!r->stats) {
        fprintf(stderr, "bed: call precompute_stats before reading\n");
        return -1;
    }
    if (b->cur_poly >= r->n_variants) return 0;  /* EOF */

    if (fread(b->raw_buf, 1, b->bytes_per_var, b->fp) != b->bytes_per_var)
        return feof(b->fp) ? 0 : -1;

    decode_raw(b->raw_buf, b->n_samples_total, out);
    if (flags) apply_flags(out, b->n_samples_total, &r->stats[b->cur_poly], flags);
    b->cur_poly++;

    /* seek past any intervening monomorphic variants */
    if (b->cur_poly < r->n_variants) {
        int next = b->poly_idx[b->cur_poly];
        int prev = b->poly_idx[b->cur_poly - 1];
        if (next != prev + 1) {
            off_t off = 3 + (off_t)next * (off_t)b->bytes_per_var;
            fseeko(b->fp, off, SEEK_SET);
        }
    }
    return 1;
}

static int bed_read_chunk(geno_reader_t *r, int k, double *out, int flags)
{
    int n = r->n_samples;
    int nread = 0;
    while (nread < k) {
        int ret = bed_read_variant(r, out + (size_t)nread * n, flags);
        if (ret <= 0) break;
        nread++;
    }
    return nread;
}

static void bed_close(geno_reader_t *r)
{
    bed_impl_t *b = r->impl;
    if (b->fp) fclose(b->fp);
    free(b->raw_buf);
    free(b->poly_idx);
    free(b);
    free(r->stats);
    free(r);
}

static const geno_ops_t bed_ops = {
    bed_precompute_stats,
    bed_read_variant,
    bed_read_chunk,
    bed_seek,
    bed_rewind,
    bed_close,
};

geno_reader_t *geno_open_bed(const char *prefix)
{
    char path[4096];

    snprintf(path, sizeof(path), "%s.fam", prefix);
    int n_samples = count_lines(path);
    if (n_samples <= 0) {
        fprintf(stderr, "bed: cannot read %s\n", path);
        return NULL;
    }

    snprintf(path, sizeof(path), "%s.bim", prefix);
    int n_variants = count_lines(path);
    if (n_variants <= 0) {
        fprintf(stderr, "bed: cannot read %s\n", path);
        return NULL;
    }

    snprintf(path, sizeof(path), "%s.bed", prefix);
    FILE *fp = fopen(path, "rb");
    if (!fp) {
        fprintf(stderr, "bed: cannot open %s: %s\n", path, strerror(errno));
        return NULL;
    }

    unsigned char magic[3];
    if (fread(magic, 1, 3, fp) != 3 || memcmp(magic, BED_MAGIC, 3) != 0) {
        fprintf(stderr, "bed: bad magic bytes in %s\n", path);
        fclose(fp);
        return NULL;
    }

    bed_impl_t *b = calloc(1, sizeof(*b));
    if (!b) { fclose(fp); return NULL; }

    b->fp               = fp;
    b->n_samples_total  = n_samples;
    b->n_variants_total = n_variants;
    b->bytes_per_var    = ((size_t)n_samples + 3) / 4;
    b->raw_buf          = malloc(b->bytes_per_var);
    if (!b->raw_buf) { fclose(fp); free(b); return NULL; }

    geno_reader_t *r = calloc(1, sizeof(*r));
    if (!r) { fclose(fp); free(b->raw_buf); free(b); return NULL; }

    r->ops        = &bed_ops;
    r->n_samples  = n_samples;
    r->n_variants = n_variants;  /* updated to poly count by precompute_stats */
    r->impl       = b;

    return r;
}
