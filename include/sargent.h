#pragma once

#include <tskit.h>
#include "geno.h"

typedef struct {
    int n;              /* number of individuals */
    int p;              /* number of markers */
    int S;              /* number of segments (local trees) */
    geno_reader_t  *gr; /* genotype reader (borrowed) */
    tsk_treeseq_t  *ts; /* tree sequence (borrowed) */
    double         *y;  /* phenotype vector [n] */
    double        *beta;   /* fixed effects [p], NULL until inferred */
    double        *ws2;    /* segment weights [S], NULL until inferred */
    double h2;          /* total heritability */
    double h2_fixed;    /* fixed-effect heritability component */
    double h2_random;   /* random-effect heritability component */
    double sigma2;      /* nongenetic (residual) variance */
} sargent_t;

sargent_t *sargent_new(geno_reader_t *gr, tsk_treeseq_t *ts, double *y);
void       sargent_free(sargent_t *sg);
