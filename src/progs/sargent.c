#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <getopt.h>
#include <tskit.h>
#include "geno.h"
#include "sargent.h"
#include "variational.h"
#include "sargent.help"
#include <phast/misc.h>

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
    int silent  = 0;
    int ts_geno = 0;

    struct option long_opts[] = {
        {"silent",  0, 0, 'z'},
        {"ts-geno", 0, 0, 'G'},
        {"version", 0, 0, 'v'},
        {"help",    0, 0, 'h'},
        {0, 0, 0, 0}
    };

    while ((c = getopt_long(argc, argv, "zGvh", long_opts, &opt_idx)) != -1) {
        switch (c) {
        case 'z': silent  = 1; break;
        case 'G': ts_geno = 1; break;
        case 'v': printf("sargent version %s\n", SARGENT_VERSION); exit(0);
        case 'h': printf("%s", HELP); exit(0);
        case '?':
            fprintf(stderr, "Bad argument.  Try 'sargent -h'.\n");
            exit(EXIT_FAILURE);
        }
    }

    int n_pos = argc - optind;
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

    FILE *pheno_fh = phast_fopen(pheno_file, "r");
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

    /* Load tree sequence (used for local-tree covariance and optionally genotypes) */
    if (!silent) fprintf(stderr, "Reading %s...\n", trees_file);
    tsk_treeseq_t ts;
    CHECK(tsk_treeseq_load(&ts, trees_file, 0));

    /* Open genotype reader */
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

    /* Precompute per-SNP statistics (allele frequencies, standardization params) */
    if (!silent) fprintf(stderr, "Precomputing variant statistics...\n");
    if (geno_precompute_stats(gr) < 0) {
        fprintf(stderr, "ERROR: precompute_stats failed.\n");
        geno_close(gr);
        tsk_treeseq_free(&ts);
        exit(EXIT_FAILURE);
    }

    if (!silent)
        fprintf(stderr, "Polymorphic variants: %d\n", gr->n_variants);

    /* Read phenotypes now that n is known */
    double *y = read_phenotypes(pheno_fh, gr->n_samples);
    fclose(pheno_fh);
    if (!y) {
        geno_close(gr);
        tsk_treeseq_free(&ts);
        exit(EXIT_FAILURE);
    }

    /* Build the sargent object and run inference */
    sargent_t *sg = sargent_new(gr, &ts, y);
    free(y);
    if (!sg) {
        fprintf(stderr, "ERROR: failed to build sargent object.\n");
        geno_close(gr);
        tsk_treeseq_free(&ts);
        exit(EXIT_FAILURE);
    }

    if (!silent) fprintf(stderr, "Running var_fixed...\n");
    var_fixed(sg);
    if (!silent) fprintf(stderr, "h2_fixed = %.4f\n", sg->h2_fixed);

    sargent_free(sg);
    geno_close(gr);
    tsk_treeseq_free(&ts);
    if (!silent) fprintf(stderr, "Done.\n");
    return 0;
}
