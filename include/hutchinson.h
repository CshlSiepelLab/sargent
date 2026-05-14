#pragma once

/* Hutchinson stochastic trace estimator for SARGENT.
 *
 * Estimates B[s] = tr(V^{-1} G_s) for all S local trees in one efficient pass:
 *
 *   1. Draw N_probe Rademacher probe vectors eps[i].
 *   2. Solve V z[i] = eps[i] for each i via CG (N_probe solves total).
 *   3. Walk the tree sequence once; for each tree s compute G_s z[i] in O(n)
 *      and accumulate eps[i]^T (G_s z[i]) into B[s].
 *
 * This avoids the naive S × N_probe CG solves. */

#include <tskit.h>
#include "conjgrad.h"

typedef struct {
    int          nprobe;       /* number of Hutchinson probes (default 50) */
    int          batch_size;   /* probes held in memory at once; 0 = all nprobe.
                                  reduce to limit RAM at the cost of extra tree passes */
    int          use_gaussian; /* 0 = Rademacher (default), 1 = standard normal */
    unsigned int seed;         /* RNG seed; 0 = use rand() state as-is */
} hutch_opts_t;

/* Estimate B[s] = tr(V^{-1} G_s) for s = 0 .. S-1.
 *
 * cg:      pre-allocated CG state (n must match).
 * Av:      matvec callback for V = sum_s ws2[s] G_s + sigma2 I.
 * av_data: opaque pointer forwarded to Av.
 * ts:      tree sequence (S trees, n diploid samples).
 * n:       number of individuals.
 * S:       number of local trees (must equal tsk_treeseq_get_num_trees(ts)).
 * B:       output array, length S (caller-allocated, overwritten).
 * opts:    tuning options; NULL for defaults (nprobe=50, Rademacher, all in one batch). */
/* Estimate B[s] = tr(V^{-1} G_s) for all S trees and, optionally,
 * B0 = tr(V^{-1}) (written to *B0_out if non-NULL).
 *
 * B0 is computed at no extra cost as (1/nprobe) sum_i eps_i^T z_i using the
 * same probe/solution pairs already needed for B[s]. */
void hutch_trace_all(cg_state_t         *cg,
                     cg_matvec_fn        Av,
                     void               *av_data,
                     const tsk_treeseq_t *ts,
                     int                 n,
                     int                 S,
                     double             *B,
                     double             *B0_out,   /* may be NULL */
                     const hutch_opts_t *opts);

/* Lower-level: given pre-computed probe matrix eps[nprobe*n] and solution
 * matrix z[nprobe*n] (row-major, each row a vector of length n), walk the
 * tree sequence and accumulate B[s] estimates.  If B0_out is non-NULL,
 * also accumulates B0 = tr(V^{-1}) = mean(eps_i^T z_i). */
void hutch_accumulate(const double        *eps,
                      const double        *z,
                      int                  nprobe,
                      int                  n,
                      const tsk_treeseq_t *ts,
                      double              *B,
                      int                  S,
                      double              *B0_out);
