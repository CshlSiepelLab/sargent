/*
 * SARGENT: Sparse ARG-based ENcoding of Traits
 *
 * Copyright (c) 2026, Cold Spring Harbor Laboratory
 * All rights reserved.
 *
 * This file is part of SARGENT and is distributed under the BSD 3-Clause License.
 * See the LICENSE file in the project root for details.
 */
#ifndef SARGENT_THREADS_H
#define SARGENT_THREADS_H

/* Establish a single thread budget across EVERY pool this binary can spin up.
 *
 * Why this exists: sargent links both libgomp (OpenMP, our own parallel loops)
 * and a BLAS (LAPACK dsyev_ in gall.c).  OpenBLAS keeps its OWN thread pool and
 * does NOT consult OMP_NUM_THREADS, so setting that variable caps our loops and
 * leaves the eigendecomposition free to seize every core.  On a shared 56-core
 * machine a single --uniform-arg run took 44 cores and drove load to 94 while
 * OMP_NUM_THREADS was set to 4.  Callers reasonably expect one knob to work.
 *
 * Resolution order (first hit wins):
 *   1. requested          (--threads N, N > 0)
 *   2. SARGENT_THREADS    (env)
 *   3. OPENBLAS_NUM_THREADS / MKL_NUM_THREADS (env; already honoured by the
 *      library, applied here too so all pools agree)
 *   4. OMP_NUM_THREADS    (env)  <-- the least-surprise fix: this now caps BLAS
 *   5. unset: leave library defaults alone (no behaviour change)
 *
 * Returns the budget applied, or 0 if defaults were left untouched.
 * Call once, early in main(), before any BLAS or OpenMP work.
 */
int sargent_init_threads(int requested, int silent);

#endif /* SARGENT_THREADS_H */
