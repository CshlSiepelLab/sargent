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
#include "threads.h"

#ifdef _OPENMP
#include <omp.h>
#endif

/* Weak references: the binary links whether or not these symbols exist, so we
 * do not take a hard dependency on any particular BLAS.  goto_set_num_threads
 * is OpenBLAS's older alias and is still exported by some builds. */
#if defined(__GNUC__) || defined(__clang__)
extern void openblas_set_num_threads(int) __attribute__((weak));
extern void goto_set_num_threads(int)     __attribute__((weak));
extern void mkl_set_num_threads(int)      __attribute__((weak));
#define SARGENT_HAVE_WEAK 1
#endif

static int env_int(const char *name)
{
    const char *v = getenv(name);
    if (!v || !*v) return 0;
    int n = atoi(v);
    return (n > 0) ? n : 0;
}

int sargent_init_threads(int requested, int silent)
{
    int n = 0;
    const char *src = NULL;

    if (requested > 0)                        { n = requested;                        src = "--threads"; }
    else if ((n = env_int("SARGENT_THREADS")))          src = "SARGENT_THREADS";
    else if ((n = env_int("OPENBLAS_NUM_THREADS")))     src = "OPENBLAS_NUM_THREADS";
    else if ((n = env_int("MKL_NUM_THREADS")))          src = "MKL_NUM_THREADS";
    else if ((n = env_int("OMP_NUM_THREADS")))          src = "OMP_NUM_THREADS";

    if (n <= 0) return 0;   /* nothing requested: leave library defaults alone */

#ifdef _OPENMP
    omp_set_num_threads(n);
#endif

#ifdef SARGENT_HAVE_WEAK
    if (openblas_set_num_threads) openblas_set_num_threads(n);
    if (goto_set_num_threads)     goto_set_num_threads(n);
    if (mkl_set_num_threads)      mkl_set_num_threads(n);
#endif

    /* Belt and braces: some BLAS builds read these only at first use, and the
     * setter may be absent from the build we linked against. */
    {
        char buf[32];
        snprintf(buf, sizeof buf, "%d", n);
        setenv("OMP_NUM_THREADS",      buf, 1);
        setenv("OPENBLAS_NUM_THREADS", buf, 1);
        setenv("GOTO_NUM_THREADS",     buf, 1);
        setenv("MKL_NUM_THREADS",      buf, 1);
    }

    if (!silent)
        fprintf(stderr, "threads: %d (from %s), applied to OpenMP and BLAS\n", n, src);
    return n;
}
