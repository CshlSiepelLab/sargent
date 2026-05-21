/*
 * SARGENT — Phase 5 microbenchmark.
 *
 * Times the two production-candidate implementations of out = G_s * in
 * across a tree sequence:
 *   - tree_cov_matvec               (current naive: per-tree two-pass)
 *   - tree_cov_matvec_incremental   (Phase 5 chunk 2+3)
 *
 * Reports wall-time per full tree-sequence walk.  Use to gauge the
 * incremental speedup before / after each chunk lands.
 *
 * Usage:
 *   bench_tree_matvec <input.trees> [n_probes=200] [seed=42]
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <tskit.h>

#include "geno.h"

static double now_seconds(void)
{
    struct timespec ts;
    clock_gettime(CLOCK_MONOTONIC, &ts);
    return ts.tv_sec + ts.tv_nsec * 1e-9;
}

/* xorshift64* + Box-Muller standard normal. */
static double rand_normal(uint64_t *state)
{
    *state ^= *state >> 12; *state ^= *state << 25; *state ^= *state >> 27;
    uint64_t r1 = *state * 0x2545F4914F6CDD1DULL;
    *state ^= *state >> 12; *state ^= *state << 25; *state ^= *state >> 27;
    uint64_t r2 = *state * 0x2545F4914F6CDD1DULL;
    double u1 = (double)(r1 >> 11) / (double)(1ULL << 53);
    double u2 = (double)(r2 >> 11) / (double)(1ULL << 53);
    if (u1 < 1e-300) u1 = 1e-300;
    return sqrt(-2.0 * log(u1)) * cos(2.0 * M_PI * u2);
}

int main(int argc, char **argv)
{
    if (argc < 2) {
        fprintf(stderr, "Usage: %s <input.trees> [n_probes=200] [seed=42]\n", argv[0]);
        return 2;
    }
    const char *trees_path = argv[1];
    int n_probes = (argc > 2) ? atoi(argv[2]) : 200;
    uint64_t seed = (argc > 3) ? (uint64_t)strtoull(argv[3], NULL, 10) : 42ULL;
    if (seed == 0) seed = 42ULL;

    tsk_treeseq_t ts;
    if (tsk_treeseq_load(&ts, trees_path, 0) < 0) {
        fprintf(stderr, "ERROR: cannot load %s\n", trees_path);
        return 2;
    }
    int *hap_to_indiv = tsk_build_hap_to_indiv(&ts);
    if (!hap_to_indiv) { tsk_treeseq_free(&ts); return 2; }

    tsk_size_t n_haps    = tsk_treeseq_get_num_samples(&ts);
    tsk_size_t n_nodes   = tsk_treeseq_get_num_nodes(&ts);
    tsk_size_t n_indiv_t = tsk_treeseq_get_num_individuals(&ts);
    int n_indiv = (n_indiv_t > 0) ? (int)n_indiv_t : (int)(n_haps / 2);
    int n_trees = (int)tsk_treeseq_get_num_trees(&ts);

    fprintf(stderr, "bench_tree_matvec: %s\n", trees_path);
    fprintf(stderr, "  n_indiv=%d  n_haps=%lld  n_nodes=%lld  n_trees=%d  n_probes=%d\n",
            n_indiv, (long long)n_haps, (long long)n_nodes, n_trees, n_probes);

    /* Scratch buffers */
    double   *in        = malloc((size_t)n_indiv * sizeof(double));
    double   *out       = malloc((size_t)n_indiv * sizeof(double));
    double   *work      = malloc(2 * n_nodes * sizeof(double));
    tsk_id_t *order_buf = malloc(n_nodes * sizeof(tsk_id_t));
    tree_cov_matvec_inc_state_t inc;
    if (!in || !out || !work || !order_buf ||
        tree_cov_matvec_inc_state_init(&inc, &ts, n_indiv, hap_to_indiv) < 0) {
        fprintf(stderr, "ERROR: out of memory\n");
        return 2;
    }

    /* Pre-generate `in` vectors so RNG cost is excluded from timing. */
    double *in_buf = malloc((size_t)n_probes * n_indiv * sizeof(double));
    for (int k = 0; k < n_probes; k++)
        for (int i = 0; i < n_indiv; i++)
            in_buf[(size_t)k * n_indiv + i] = rand_normal(&seed);

    /* sink: accumulate a hash so the compiler can't drop the work. */
    double sink_naive = 0.0, sink_inc = 0.0;

    /* --- Time naive --- */
    double t0 = now_seconds();
    for (int k = 0; k < n_probes; k++) {
        const double *in_k = in_buf + (size_t)k * n_indiv;
        tsk_tree_t tree;
        if (tsk_tree_init(&tree, &ts, 0) < 0) break;
        for (int ret = tsk_tree_first(&tree); ret == TSK_TREE_OK;
             ret = tsk_tree_next(&tree)) {
            tree_cov_matvec(out, in_k, n_indiv, &tree, work, order_buf, hap_to_indiv);
            for (int i = 0; i < n_indiv; i++) sink_naive += out[i];
        }
        tsk_tree_free(&tree);
    }
    double t_naive = now_seconds() - t0;

    /* --- Time incremental --- */
    t0 = now_seconds();
    for (int k = 0; k < n_probes; k++) {
        const double *in_k = in_buf + (size_t)k * n_indiv;
        tree_cov_matvec_inc_state_reset(&inc, in_k);
        tsk_tree_t tree;
        if (tsk_tree_init(&tree, &ts, 0) < 0) break;
        for (int ret = tsk_tree_first(&tree); ret == TSK_TREE_OK;
             ret = tsk_tree_next(&tree)) {
            tree_cov_matvec_incremental(out, &inc, &tree);
            for (int i = 0; i < n_indiv; i++) sink_inc += out[i];
        }
        tsk_tree_free(&tree);
    }
    double t_inc = now_seconds() - t0;

    fprintf(stderr, "\n  naive       : %.3f s total  ->  %.3f ms / sequence walk\n",
            t_naive, 1000.0 * t_naive / n_probes);
    fprintf(stderr, "  incremental : %.3f s total  ->  %.3f ms / sequence walk\n",
            t_inc,   1000.0 * t_inc   / n_probes);
    fprintf(stderr, "  speedup     : %.2fx\n", t_naive / t_inc);
    fprintf(stderr, "  sink check  : naive=%.4e  inc=%.4e  diff=%.2e\n",
            sink_naive, sink_inc, fabs(sink_naive - sink_inc));

    free(in); free(out); free(work); free(order_buf); free(in_buf);
    tree_cov_matvec_inc_state_free(&inc);
    free(hap_to_indiv);
    tsk_treeseq_free(&ts);
    return 0;
}
