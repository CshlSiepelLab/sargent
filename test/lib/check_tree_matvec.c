/*
 * SARGENT — Phase 5 paranoia harness for tree-cov matvec implementations.
 *
 * Loads a tree sequence, generates K random Gaussian probe vectors, and for
 * each (tree, probe) compares two independent implementations of
 *   out = G_s * in
 * (where G_s[i,j] = 0.5 * (time_root - time_MRCA(haps_of_i, haps_of_j))):
 *   1. tree_cov_matvec           (production: post-/pre-order two-pass)
 *   2. tree_cov_matvec_pairwise  (reference:  direct pairwise MRCA sum)
 *
 * Asserts L_inf agreement < tol.  Reports max abs / max rel diff across
 * all (tree, probe) pairs.  Exits 0 on success, 1 on first divergence.
 *
 * Usage:
 *   check_tree_matvec <input.trees> [n_probes] [max_trees] [tol] [seed]
 *
 * tol is interpreted as a MIXED tolerance: failure when
 *     diff > tol * (1 + max(|naive|, |ref|))
 * which is robust across magnitudes (branch lengths can run 10^6).
 *
 * Defaults: n_probes=100, max_trees=0 (all), tol=1e-9, seed=42.
 (1e-9 absorbs the worst-case O(n_terms * ulp) accumulation across the
 two algorithms; in practice we observe ~1e-10 on debug cases.)
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <tskit.h>

#include "geno.h"   /* tree_cov_matvec, tree_cov_matvec_pairwise, tsk_build_hap_to_indiv */

/* Box-Muller standard normal, using a 64-bit xorshift PRNG. */
static double rand_normal(uint64_t *state)
{
    /* xorshift64* */
    *state ^= *state >> 12;
    *state ^= *state << 25;
    *state ^= *state >> 27;
    uint64_t r1 = *state * 0x2545F4914F6CDD1DULL;

    *state ^= *state >> 12;
    *state ^= *state << 25;
    *state ^= *state >> 27;
    uint64_t r2 = *state * 0x2545F4914F6CDD1DULL;

    /* Map two uint64 to two uniforms in (0, 1] then Box-Muller. */
    double u1 = (double)(r1 >> 11) / (double)(1ULL << 53);
    double u2 = (double)(r2 >> 11) / (double)(1ULL << 53);
    if (u1 < 1e-300) u1 = 1e-300;
    return sqrt(-2.0 * log(u1)) * cos(2.0 * M_PI * u2);
}

int main(int argc, char **argv)
{
    if (argc < 2) {
        fprintf(stderr,
                "Usage: %s <input.trees> [n_probes=100] [max_trees=0] [tol=1e-10] [seed=42]\n",
                argv[0]);
        return 2;
    }
    const char *trees_path = argv[1];
    int    n_probes  = (argc > 2) ? atoi(argv[2]) : 100;
    int    max_trees = (argc > 3) ? atoi(argv[3]) : 0;
    double tol       = (argc > 4) ? atof(argv[4]) : 1e-9;
    uint64_t seed    = (argc > 5) ? (uint64_t)strtoull(argv[5], NULL, 10) : 42ULL;
    if (seed == 0) seed = 42ULL;

    tsk_treeseq_t ts;
    int ret = tsk_treeseq_load(&ts, trees_path, 0);
    if (ret < 0) {
        fprintf(stderr, "ERROR: cannot load %s: %s\n", trees_path, tsk_strerror(ret));
        return 2;
    }

    int *hap_to_indiv = tsk_build_hap_to_indiv(&ts);
    if (!hap_to_indiv) {
        fprintf(stderr, "ERROR: tsk_build_hap_to_indiv failed\n");
        tsk_treeseq_free(&ts);
        return 2;
    }

    tsk_size_t n_haps    = tsk_treeseq_get_num_samples(&ts);
    tsk_size_t n_nodes   = tsk_treeseq_get_num_nodes(&ts);
    tsk_size_t n_indiv_t = tsk_treeseq_get_num_individuals(&ts);
    int n_indiv = (n_indiv_t > 0) ? (int)n_indiv_t : (int)(n_haps / 2);

    int S_total = (int)tsk_treeseq_get_num_trees(&ts);
    int S       = (max_trees > 0 && max_trees < S_total) ? max_trees : S_total;

    fprintf(stderr,
            "check_tree_matvec: %s\n"
            "  n_indiv=%d  n_haps=%lld  n_nodes=%lld  n_trees=%d  (checking first %d)\n"
            "  n_probes=%d  tol=%g  seed=%llu\n",
            trees_path, n_indiv, (long long)n_haps, (long long)n_nodes,
            S_total, S, n_probes, tol, (unsigned long long)seed);

    /* Scratch buffers */
    double   *in        = malloc((size_t)n_indiv * sizeof(double));
    double   *out_naive = malloc((size_t)n_indiv * sizeof(double));
    double   *out_ref   = malloc((size_t)n_indiv * sizeof(double));
    double   *work      = malloc(2 * n_nodes * sizeof(double));
    tsk_id_t *order_buf = malloc(n_nodes * sizeof(tsk_id_t));
    int      *mark      = malloc(n_nodes * sizeof(int));

    if (!in || !out_naive || !out_ref || !work || !order_buf || !mark) {
        fprintf(stderr, "ERROR: out of memory\n");
        return 2;
    }

    double max_abs_diff = 0.0;
    double max_rel_diff = 0.0;
    int    n_compare    = 0;
    int    n_fail       = 0;
    int    error        = 0;

    for (int k = 0; k < n_probes; k++) {
        /* Fresh probe each iteration */
        for (int i = 0; i < n_indiv; i++) in[i] = rand_normal(&seed);

        tsk_tree_t tree;
        ret = tsk_tree_init(&tree, &ts, 0);
        if (ret < 0) {
            fprintf(stderr, "ERROR: tsk_tree_init: %s\n", tsk_strerror(ret));
            error = 1; break;
        }
        int s = 0;
        for (ret = tsk_tree_first(&tree); ret == TSK_TREE_OK && s < S;
             ret = tsk_tree_next(&tree), s++) {

            tree_cov_matvec(out_naive, in, n_indiv, &tree, work, order_buf, hap_to_indiv);
            tree_cov_matvec_pairwise(out_ref, in, n_indiv, &tree, mark, hap_to_indiv);

            for (int i = 0; i < n_indiv; i++) {
                double a = out_naive[i];
                double b = out_ref[i];
                double d = fabs(a - b);
                double mag = fmax(fabs(a), fabs(b));
                if (d > max_abs_diff) max_abs_diff = d;
                if (mag > 0.0) {
                    double r = d / mag;
                    if (r > max_rel_diff) max_rel_diff = r;
                }
                /* Mixed tolerance: tol * (1 + mag).  Don't bail -- just count. */
                if (d > tol * (1.0 + mag)) {
                    if (n_fail < 3) {
                        fprintf(stderr,
                                "  diff probe=%d tree=%d indiv=%d  naive=%.17g  ref=%.17g  "
                                "abs=%.3e  rel=%.3e\n",
                                k, s, i, a, b, d, (mag > 0.0) ? d / mag : 0.0);
                    }
                    n_fail++;
                }
            }
            n_compare++;
        }
        if (ret < 0) {
            fprintf(stderr, "ERROR: tree iteration: %s\n", tsk_strerror(ret));
            error = 1;
        }
        tsk_tree_free(&tree);
        if (error) break;
    }

    int fail = error || (n_fail > 0);
    fprintf(stderr,
            "  n_compare=%d  n_fail=%d  max_abs_diff=%.3e  max_rel_diff=%.3e  tol=%.1e  -> %s\n",
            n_compare, n_fail, max_abs_diff, max_rel_diff, tol, fail ? "FAIL" : "OK");

    free(in); free(out_naive); free(out_ref);
    free(work); free(order_buf); free(mark);
    free(hap_to_indiv);
    tsk_treeseq_free(&ts);
    return fail ? 1 : 0;
}
