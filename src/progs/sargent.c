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
#include <getopt.h>
#include <tskit.h>
#include "geno.h"
#include "sargent.h"
#include "variational.h"
#include "sargent.help.h"

/* Read n phenotype values from f.  Each non-comment line may have one or more
   whitespace-separated fields; the last field is taken as the phenotype. */
static double *read_phenotypes(FILE *f, int n)
{
    double *y = malloc((size_t)n * sizeof(double));
    if (!y) return NULL;

    char line[4096];
    int i = 0;
    while (i < n && fgets(line, sizeof(line), f)) {
        if (line[0] == '#') continue;
        char *last = NULL, *tok = strtok(line, " \t\n\r");
        while (tok) { last = tok; tok = strtok(NULL, " \t\n\r"); }
        if (!last) continue;
        char *end;
        y[i] = strtod(last, &end);
        if (end == last) {
            fprintf(stderr, "ERROR: bad phenotype value on line %d\n", i + 1);
            free(y); return NULL;
        }
        i++;
    }
    if (i < n) {
        fprintf(stderr, "ERROR: expected %d phenotypes, got %d\n", n, i);
        free(y); return NULL;
    }
    return y;
}

static void write_beta_table(sargent_t *sg, const char *outroot, int silent)
{
    int p = sg->p;
    int K = sg->K;

    char fname[4096];
    snprintf(fname, sizeof(fname), "%s.beta.tsv", outroot);
    FILE *f = fopen(fname, "w");
    if (!f) {
        fprintf(stderr, "ERROR: cannot open %s for writing\n", fname);
        return;
    }

    double *gamma_all = var_fixed_posteriors(sg);
    if (!gamma_all) {
        fprintf(stderr, "WARNING: could not compute posteriors; skipping %s\n", fname);
        fclose(f);
        return;
    }

    for (int k = 0; k < K; k++)
        fprintf(f, "gamma_%d\t", k + 1);
    fprintf(f, "E_beta\n");

    for (int j = 0; j < p; j++) {
        for (int k = 0; k < K; k++)
            fprintf(f, "%.6g\t", gamma_all[j * K + k]);
        fprintf(f, "%.6g\n", sg->beta[j]);
    }

    free(gamma_all);
    fclose(f);
    if (!silent) fprintf(stderr, "Wrote %s\n", fname);
}

static void write_w2_table(sargent_t *sg, const char *outroot, int silent)
{
    int S = sg->S;
    int M = sg->M;

    char fname[4096];
    snprintf(fname, sizeof(fname), "%s.w2.tsv", outroot);
    FILE *f = fopen(fname, "w");
    if (!f) {
        fprintf(stderr, "ERROR: cannot open %s for writing\n", fname);
        return;
    }

    for (int m = 0; m < M; m++)
        fprintf(f, "rho_%d\t", m + 1);
    fprintf(f, "E_ws2\n");

    for (int s = 0; s < S; s++) {
        for (int m = 0; m < M; m++)
            fprintf(f, "%.6g\t", sg->rho[s * M + m]);
        fprintf(f, "%.6g\n", sg->ws2[s]);
    }

    fclose(f);
    if (!silent) fprintf(stderr, "Wrote %s\n", fname);
}

#define CHECK(ret) do {                                                 \
    if ((ret) < 0) {                                                    \
        fprintf(stderr, "tskit error at %s:%d: %s\n",                  \
                __FILE__, __LINE__, tsk_strerror(ret));                 \
        exit(EXIT_FAILURE);                                             \
    }                                                                   \
} while (0)

int main(int argc, char *argv[])
{
    signed char c;
    int opt_idx;
    int silent      = 0;
    int ts_geno     = 0;
    int fixed_only  = 0;
    int random_only = 0;
    int horseshoe   = 0;
    int nprobe      = 0;   /* 0 = use sargent_t default (10) */
    double pi0      = -1.0; /* < 0 means unset */
    const char *outroot = "sargent";

    struct option long_opts[] = {
        {"silent",      0, 0, 'z'},
        {"ts-geno",     0, 0, 'G'},
        {"fixed-only",  0, 0, 'F'},
        {"random-only", 0, 0, 'R'},
        {"horseshoe",   0, 0, 'H'},
        {"pi0",         1, 0, 'P'},
        {"outroot",     1, 0, 'o'},
        {"nprobe",      1, 0, 'N'},
        {"version",     0, 0, 'v'},
        {"help",        0, 0, 'h'},
        {0, 0, 0, 0}
    };

    while ((c = getopt_long(argc, argv, "zGFRHo:N:P:vh", long_opts, &opt_idx)) != -1) {
        switch (c) {
        case 'z': silent      = 1;                break;
        case 'G': ts_geno     = 1;                break;
        case 'F': fixed_only  = 1;                break;
        case 'R': random_only = 1;                break;
        case 'H': horseshoe   = 1;                break;
        case 'P': pi0         = atof(optarg);     break;
        case 'o': outroot     = optarg;           break;
        case 'N': nprobe      = atoi(optarg);     break;
        case 'v': printf("sargent version %s\n", SARGENT_VERSION); exit(0);
        case 'h': printf("%s", HELP); exit(0);
        case '?':
            fprintf(stderr, "Bad argument.  Try 'sargent -h'.\n");
            exit(EXIT_FAILURE);
        }
    }

    if (fixed_only && random_only) {
        fprintf(stderr, "ERROR: --fixed-only and --random-only are mutually exclusive.\n");
        exit(EXIT_FAILURE);
    }

    int n_pos   = argc - optind;
    int need_pos = ts_geno ? 2 : 3;
    if (n_pos != need_pos) {
        if (ts_geno)
            fprintf(stderr,
                    "ERROR: expected <phenotypes> <input.trees> with --ts-geno.\n"
                    "Try 'sargent -h'.\n");
        else
            fprintf(stderr,
                    "ERROR: expected <phenotypes> <input.trees> <geno_prefix>.\n"
                    "Use --ts-geno to obtain genotypes from the tree sequence.\n"
                    "Try 'sargent -h'.\n");
        exit(EXIT_FAILURE);
    }

    const char *pheno_file  = argv[optind];
    const char *trees_file  = argv[optind + 1];
    const char *geno_prefix = ts_geno ? NULL : argv[optind + 2];

    FILE *pheno_fh = fopen(pheno_file, "r");
    if (!pheno_fh) {
        fprintf(stderr, "ERROR: cannot open %s\n", pheno_file);
        exit(EXIT_FAILURE);
    }

    if (!silent) {
        fprintf(stderr, "Phenotypes  : %s\n", pheno_file);
        fprintf(stderr, "Tree seq    : %s\n", trees_file);
        if (ts_geno)
            fprintf(stderr, "Genotypes   : from tree sequence\n");
        else
            fprintf(stderr, "Genotypes   : %s.bed\n", geno_prefix);
    }

    if (!silent) fprintf(stderr, "Reading %s...\n", trees_file);
    tsk_treeseq_t ts;
    CHECK(tsk_treeseq_load(&ts, trees_file, 0));

    geno_reader_t *gr = ts_geno ? geno_open_trees(&ts)
                                : geno_open_bed(geno_prefix);
    if (!gr) {
        fprintf(stderr, "ERROR: failed to open genotype source.\n");
        tsk_treeseq_free(&ts);
        exit(EXIT_FAILURE);
    }

    if (!silent)
        fprintf(stderr, "Samples: %d  Variants (pre-filter): %d\n",
                gr->n_samples, gr->n_variants);

    if (!silent) fprintf(stderr, "Precomputing variant statistics...\n");
    if (geno_precompute_stats(gr) < 0) {
        fprintf(stderr, "ERROR: precompute_stats failed.\n");
        geno_close(gr);
        tsk_treeseq_free(&ts);
        exit(EXIT_FAILURE);
    }

    if (!silent)
        fprintf(stderr, "Polymorphic variants: %d\n", gr->n_variants);

    double *y = read_phenotypes(pheno_fh, gr->n_samples);
    fclose(pheno_fh);
    if (!y) {
        geno_close(gr);
        tsk_treeseq_free(&ts);
        exit(EXIT_FAILURE);
    }

    sargent_t *sg = sargent_new(gr, &ts, y);
    free(y);
    if (!sg) {
        fprintf(stderr, "ERROR: failed to build sargent object.\n");
        geno_close(gr);
        tsk_treeseq_free(&ts);
        exit(EXIT_FAILURE);
    }
    if (nprobe > 0)  sg->nprobe = nprobe;
    if (horseshoe)   sg->prior_type = 1;
    if (pi0 > 0.0 && pi0 < 1.0 && !horseshoe) {
        /* Beta(1, pi0/(1-pi0)) hyperprior: prior mean of phi = 1-pi0.
         * E.g. pi0=0.9 → Beta(1,9), prior mean phi=0.1 (10% active). */
        sg->phi_alpha = 1.0;
        sg->phi_beta  = pi0 / (1.0 - pi0);
    }

    /* ---- Stage 1: fixed effects ---- */
    if (!random_only) {
        if (!silent)
            fprintf(stderr, "Stage 1: inferring fixed effects (%d SNPs, K=%d)...\n",
                    sg->p, sg->K);
        var_fixed(sg, silent);
        if (!silent)
            fprintf(stderr, "Stage 1 done: h2_fixed=%.4f  sigma2=%.4f\n",
                    sg->h2_fixed, sg->sigma2);
        write_beta_table(sg, outroot, silent);
    }

    /* ---- Stage 2: random effects ---- */
    if (!fixed_only) {
        if (!silent)
            fprintf(stderr, "Stage 2: inferring segment weights (%d trees, prior=%s)...\n",
                    sg->S, sg->prior_type == 1 ? "horseshoe" : "beta-bernoulli");
        var_random(sg, silent);
        if (!silent)
            fprintf(stderr, "Stage 2 done: h2_random=%.4f  sigma2=%.4f\n",
                    sg->h2_random, sg->sigma2);
        write_w2_table(sg, outroot, silent);
    }

    /* ---- Summary statistics ---- */
    if (!silent) {
        fprintf(stderr, "\nSummary:\n");
        fprintf(stderr, "  h2        : %.4f\n", sg->h2);
        fprintf(stderr, "  h2.fixed  : %.4f\n", sg->h2_fixed);
        fprintf(stderr, "  h2.random : %.4f\n", sg->h2_random);
        fprintf(stderr, "  sigma2    : %.4f\n", sg->sigma2);
    }

    sargent_free(sg);
    geno_close(gr);
    tsk_treeseq_free(&ts);
    if (!silent) fprintf(stderr, "Done.\n");
    return 0;
}
