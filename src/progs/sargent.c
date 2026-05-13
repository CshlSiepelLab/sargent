#include <stdlib.h>
#include <stdio.h>
#include <getopt.h>
#include <tskit.h>
#include "geno.h"
#include "sargent.help"
#include <phast/misc.h>

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
    fclose(pheno_fh);

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

    /* placeholder: analysis goes here */

    geno_close(gr);
    tsk_treeseq_free(&ts);
    if (!silent) fprintf(stderr, "Done.\n");
    return 0;
}
