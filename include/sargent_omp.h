/*
 * SARGENT: Sparse ARG-based ENcoding of Traits
 *
 * Copyright (c) 2026, Cold Spring Harbor Laboratory
 * All rights reserved.
 *
 * This file is part of SARGENT and is distributed under the BSD 3-Clause License.
 * See the LICENSE file in the project root for details.
 */

#pragma once

/*
 * OpenMP compatibility shim.
 *
 * The CMake build sets SARGENT_HAS_OPENMP when find_package(OpenMP) succeeds
 * AND the target is linked against OpenMP::OpenMP_C; the compiler then also
 * defines _OPENMP whenever the OpenMP runtime is in use.  We guard
 * <omp.h> on the standard _OPENMP macro and fall through to no-op accessors
 * if OpenMP is not available, so library code can call
 * SARGENT_OMP_GET_THREAD_NUM() / SARGENT_OMP_GET_MAX_THREADS() unconditionally.
 *
 * #pragma omp directives are silently ignored by compilers without
 * OpenMP support (the standard says: "unknown #pragma" → warning at most),
 * so no #ifdef is needed around them.
 */

#ifdef _OPENMP
  #include <omp.h>
  #define SARGENT_OMP_GET_MAX_THREADS()  omp_get_max_threads()
  #define SARGENT_OMP_GET_THREAD_NUM()   omp_get_thread_num()
  #define SARGENT_OMP_IN_PARALLEL()      omp_in_parallel()
#else
  #define SARGENT_OMP_GET_MAX_THREADS()  1
  #define SARGENT_OMP_GET_THREAD_NUM()   0
  #define SARGENT_OMP_IN_PARALLEL()      0
#endif
