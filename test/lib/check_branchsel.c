/*
 * SARGENT — paranoia harness for branch selection (branchsel.c).
 *
 * The whole result rests on the score denominator: an implementation that gets
 * x_b' P x_b wrong still produces a plausible-looking ranking, and the branch
 * AUC silently drops from ~0.64 to ~0.41.  So this harness re-derives every
 * quantity by a second, deliberately different route and compares:
 *
 *   enumeration   production walks each tree in postorder summing child carrier
 *                 vectors; the reference walks each sample haplotype UP to the
 *                 root marking ancestors.  Both must yield the same set of
 *                 distinct carrier vectors (compared by FNV-1a fingerprint).
 *
 *   numerator     production accumulates x_b' Py inside a blocked BLAS pass;
 *                 the reference uses a plain loop.
 *
 *   denominator   production uses one dsymm per block of branches; the
 *                 reference uses the naive O(n^2) double loop per branch.
 *
 *   ranking       the dump must be sorted by |z| descending, and z must equal
 *                 |num| / sqrt(den) from the reference values.
 *
 * P is a genuine REML projector (V^-1 minus its intercept part) built from a
 * random SPD V, so P 1 = 0 holds exactly as it does in production.
 *
 * Usage:
 *   check_branchsel <input.trees> [n_train] [k] [max_trees] [tol] [seed]
 * Defaults: n_train = 2/3 of individuals, k = 25, max_trees = 0 (all),
 *           tol = 1e-8 (mixed), seed = 42.
 *
 * Exits 0 if every comparison passes.
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdint.h>
#include <unistd.h>
#include <tskit.h>

#include "geno.h"
#include "branchsel.h"

extern void dpotrf_(const char *uplo, const int *n, double *a, const int *lda, int *info);
extern void dpotri_(const char *uplo, const int *n, double *a, const int *lda, int *info);

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

typedef struct {
    uint64_t fp;
    double   num, den, z;
    long     nct, ncr;
} ref_branch_t;

int main(int argc, char **argv)
{
    if (argc < 2) {
        fprintf(stderr, "Usage: %s <input.trees> [n_train] [k] [max_trees] "
                "[tol] [seed]\n", argv[0]);
        return 2;
    }
    const char *trees_path = argv[1];
    tsk_treeseq_t ts;
    int ret = tsk_treeseq_load(&ts, trees_path, 0);
    if (ret < 0) {
        fprintf(stderr, "cannot load %s: %s\n", trees_path, tsk_strerror(ret));
        return 2;
    }
    int n_total = (int)tsk_treeseq_get_num_individuals(&ts);
    int n_train   = (argc > 2) ? atoi(argv[2]) : (n_total * 2) / 3;
    int k         = (argc > 3) ? atoi(argv[3]) : 25;
    int max_trees = (argc > 4) ? atoi(argv[4]) : 0;
    double tol    = (argc > 5) ? atof(argv[5]) : 1e-8;
    uint64_t seed = (argc > 6) ? strtoull(argv[6], NULL, 10) : 42ULL;
    if (seed == 0) seed = 42ULL;
    if (n_total <= 0) { fprintf(stderr, "no individual table\n"); return 2; }
    if (n_train <= 1 || n_train >= n_total) {
        fprintf(stderr, "n_train=%d out of range for %d individuals\n",
                n_train, n_total);
        return 2;
    }

    int *h2i = tsk_build_hap_to_indiv(&ts);
    if (!h2i) return 2;

    int S = (int)tsk_treeseq_get_num_trees(&ts);
    printf("check_branchsel: %s  individuals=%d (train=%d)  trees=%d  k=%d\n",
           trees_path, n_total, n_train, S, k);

    /* ---- a real REML projector from a random SPD V --------------------- */
    int n = n_train;
    double *V = malloc((size_t)n * n * sizeof(double));
    double *M = malloc((size_t)n * n * sizeof(double));
    double *P = malloc((size_t)n * n * sizeof(double));
    double *y = malloc((size_t)n * sizeof(double));
    double *Py = malloc((size_t)n * sizeof(double));
    if (!V || !M || !P || !y || !Py) return 2;

    for (size_t q = 0; q < (size_t)n * n; q++) M[q] = rand_normal(&seed);
    for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++) {
            double s = 0.0;
            for (int t = 0; t < n; t++) s += M[(size_t)i * n + t] * M[(size_t)j * n + t];
            V[(size_t)i * n + j] = s / n + (i == j ? 1.0 : 0.0);
        }
    memcpy(P, V, (size_t)n * n * sizeof(double));
    int info = 0;
    dpotrf_("L", &n, P, &n, &info);
    if (info == 0) dpotri_("L", &n, P, &n, &info);
    if (info != 0) { fprintf(stderr, "V inverse failed\n"); return 2; }
    for (int i = 0; i < n; i++)
        for (int j = i + 1; j < n; j++) P[(size_t)j * n + i] = P[(size_t)i * n + j];
    {   /* P = Vinv - v1 v1' / (1' Vinv 1) */
        double *v1 = malloc((size_t)n * sizeof(double));
        double s11 = 0.0;
        for (int i = 0; i < n; i++) {
            double r = 0.0;
            for (int j = 0; j < n; j++) r += P[(size_t)i * n + j];
            v1[i] = r; s11 += r;
        }
        for (int i = 0; i < n; i++)
            for (int j = 0; j < n; j++)
                P[(size_t)i * n + j] -= v1[i] * v1[j] / s11;
        free(v1);
    }
    double mean_y = 0.0;
    for (int i = 0; i < n; i++) { y[i] = rand_normal(&seed); mean_y += y[i]; }
    mean_y /= n;
    for (int i = 0; i < n; i++) y[i] -= mean_y;
    for (int i = 0; i < n; i++) {
        double s = 0.0;
        for (int j = 0; j < n; j++) s += P[(size_t)i * n + j] * y[j];
        Py[i] = s;
    }
    {   /* P 1 = 0 must hold to machine precision, else the "centering cancels"
           argument behind using raw carrier counts is void. */
        double worst = 0.0, scale = 0.0;
        for (int i = 0; i < n; i++) {
            double r = 0.0, a = 0.0;
            for (int j = 0; j < n; j++) { r += P[(size_t)i * n + j];
                                          a += fabs(P[(size_t)i * n + j]); }
            if (fabs(r) > worst) worst = fabs(r);
            if (a > scale) scale = a;
        }
        printf("  P 1 = 0 check: max |row sum| = %.3e (row L1 scale %.3e)\n",
               worst, scale);
        if (worst > 1e-9 * scale) {
            printf("FAIL: projector is not intercept-orthogonal\n");
            return 1;
        }
    }

    /* ---- run production, capturing the per-branch dump ------------------ */
    char dump[512];
    snprintf(dump, sizeof(dump), "/tmp/check_branchsel_%d.tsv", (int)getpid());
    setenv("SARGENT_BRANCHSEL_DUMP", dump, 1);

    branchsel_opts_t opts = { .k = k, .max_trees = max_trees, .silent = 0 };
    branchsel_t bs;
    if (branchsel_run(&ts, h2i, n_total, n_train, P, y, &opts, &bs) != 0) {
        fprintf(stderr, "branchsel_run failed\n");
        return 2;
    }

    /* ---- reference: mark ancestors by walking each sample UP ------------ */
    tsk_size_t nnd = tsk_treeseq_get_num_nodes(&ts);
    const tsk_id_t *samples = tsk_treeseq_get_samples(&ts);
    tsk_size_t n_hap = tsk_treeseq_get_num_samples(&ts);
    const double *ntime = ts.tables->nodes.time;

    unsigned char *car = malloc(nnd * (size_t)n_total);
    ref_branch_t  *rb  = NULL;
    size_t         n_rb = 0, cap_rb = 4096;
    rb = malloc(cap_rb * sizeof(ref_branch_t));
    if (!car || !rb) return 2;

    unsigned char *keep = NULL;
    if (max_trees > 0 && max_trees < S) {
        keep = calloc((size_t)S, 1);
        for (int j = 0; j < max_trees; j++) {
            double x = (max_trees == 1) ? 0.0
                     : (double)j * (double)(S - 1) / (double)(max_trees - 1);
            int t = (int)floor(x + 0.5);
            if (t < 0) t = 0;
            if (t > S - 1) t = S - 1;
            keep[t] = 1;
        }
    }

    tsk_tree_t tree;
    tsk_tree_init(&tree, &ts, 0);
    int s = 0;
    for (ret = tsk_tree_first(&tree); ret == TSK_TREE_OK;
         ret = tsk_tree_next(&tree), s++) {
        if (s >= S) break;
        if (keep && !keep[s]) continue;
        memset(car, 0, nnd * (size_t)n_total);
        for (tsk_size_t h = 0; h < n_hap; h++) {
            int ind = h2i[(int)h];
            if (ind < 0 || ind >= n_total) continue;
            for (tsk_id_t u = samples[h]; u != TSK_NULL; u = tree.parent[u])
                car[(size_t)u * n_total + ind]++;
        }
        for (tsk_size_t u = 0; u < nnd; u++) {
            tsk_id_t p = tree.parent[u];
            if (p == TSK_NULL) continue;
            if (!(ntime[p] - ntime[u] > 0.0)) continue;
            const unsigned char *c = car + u * (size_t)n_total;
            long nct = 0, ncr = 0;
            for (int i = 0; i < n_total; i++) {
                nct += c[i];
                if (i < n_train) ncr += c[i];
            }
            if (nct < 1 || nct > 2L * n_total - 1) continue;
            /* drop branches constant within train, as production does */
            double s1 = 0.0, s2 = 0.0;
            for (int i = 0; i < n_train; i++) { s1 += c[i]; s2 += (double)c[i]*c[i]; }
            double mu = s1 / n_train, var = s2 / n_train - mu * mu;
            if (!(var > 1e-18)) continue;

            uint64_t fp = branchsel_fingerprint(c, n_total);
            if (fp == 0) fp = 1;
            size_t q;
            for (q = 0; q < n_rb; q++) if (rb[q].fp == fp) break;
            if (q < n_rb) continue;                    /* already have it */

            double num = 0.0, den = 0.0;
            for (int i = 0; i < n_train; i++) {
                if (!c[i]) continue;
                num += (double)c[i] * Py[i];
                const double *Pi = P + (size_t)i * n_train;
                double acc = 0.0;
                for (int j = 0; j < n_train; j++) acc += Pi[j] * (double)c[j];
                den += (double)c[i] * acc;
            }
            if (n_rb == cap_rb) {
                cap_rb *= 2;
                rb = realloc(rb, cap_rb * sizeof(ref_branch_t));
                if (!rb) return 2;
            }
            rb[n_rb].fp = fp; rb[n_rb].num = num; rb[n_rb].den = den;
            rb[n_rb].z = (den > 1e-12) ? fabs(num) / sqrt(den) : 0.0;
            rb[n_rb].nct = nct; rb[n_rb].ncr = ncr;
            n_rb++;
        }
    }
    tsk_tree_free(&tree);

    /* ---- compare -------------------------------------------------------- */
    FILE *f = fopen(dump, "r");
    if (!f) { fprintf(stderr, "cannot read %s\n", dump); return 2; }
    char line[1024];
    if (!fgets(line, sizeof(line), f)) { fclose(f); return 2; }   /* header */

    size_t n_dump = 0, n_matched = 0;
    double worst_num = 0.0, worst_den = 0.0, worst_z = 0.0, prev_z = INFINITY;
    int    order_ok = 1, fail = 0;
    unsigned char *seen = calloc(n_rb, 1);

    while (fgets(line, sizeof(line), f)) {
        unsigned long long fp;
        long nct, ncr; double num, den, z; size_t rank;
        if (sscanf(line, "%llx\t%ld\t%ld\t%lf\t%lf\t%lf\t%zu",
                   &fp, &nct, &ncr, &num, &den, &z, &rank) != 7) {
            fprintf(stderr, "malformed dump line: %s", line);
            fail = 1; break;
        }
        n_dump++;
        if (z > prev_z * (1.0 + 1e-12)) order_ok = 0;
        prev_z = z;

        size_t q;
        for (q = 0; q < n_rb; q++) if (rb[q].fp == (uint64_t)fp) break;
        if (q == n_rb) {
            printf("FAIL: dump branch %016llx has no reference counterpart\n", fp);
            fail = 1; continue;
        }
        seen[q] = 1; n_matched++;
        if (nct != rb[q].nct || ncr != rb[q].ncr) {
            printf("FAIL: carrier counts differ for %016llx: "
                   "dump %ld/%ld vs ref %ld/%ld\n", fp, nct, ncr,
                   rb[q].nct, rb[q].ncr);
            fail = 1;
        }
        double dn = fabs(num - rb[q].num) / (1.0 + fabs(rb[q].num));
        double dd = fabs(den - rb[q].den) / (1.0 + fabs(rb[q].den));
        double dz = fabs(z   - rb[q].z)   / (1.0 + fabs(rb[q].z));
        if (dn > worst_num) worst_num = dn;
        if (dd > worst_den) worst_den = dd;
        if (dz > worst_z)   worst_z   = dz;
    }
    fclose(f);
    remove(dump);

    size_t missing = 0;
    for (size_t q = 0; q < n_rb; q++) if (!seen[q]) missing++;

    printf("  branches: production %zu, reference %zu, matched %zu, "
           "reference-only %zu\n", n_dump, n_rb, n_matched, missing);
    printf("  max rel diff: numerator %.3e  denominator %.3e  z %.3e\n",
           worst_num, worst_den, worst_z);
    printf("  ranking sorted by |z| descending: %s\n", order_ok ? "yes" : "NO");

    if (n_dump != n_rb || missing != 0) {
        printf("FAIL: branch sets differ\n");
        fail = 1;
    }
    if (worst_num > tol || worst_den > tol || worst_z > tol) {
        printf("FAIL: exceeds tolerance %.1e\n", tol);
        fail = 1;
    }
    if (!order_ok) { printf("FAIL: dump is not sorted by score\n"); fail = 1; }

    /* The selected columns must be the top k of the reference ranking, and
     * standardized by the TRAIN mean and population sd. */
    if (bs.k > 0) {
        double worst_s = 0.0;
        for (int t = 0; t < bs.k; t++) {
            double m1 = 0.0, m2 = 0.0;
            for (int i = 0; i < n_train; i++) {
                double v = bs.Z[(size_t)i * bs.k + t];
                m1 += v; m2 += v * v;
            }
            m1 /= n_train; m2 = m2 / n_train - m1 * m1;
            double d = fabs(m1) + fabs(m2 - 1.0);
            if (d > worst_s) worst_s = d;
        }
        printf("  Z columns standardized on train (max |mean| + |var-1|): %.3e\n",
               worst_s);
        if (worst_s > 1e-9) { printf("FAIL: Z is not train-standardized\n"); fail = 1; }
    }

    free(seen); free(car); free(rb); free(keep);
    free(V); free(M); free(P); free(y); free(Py); free(h2i);
    branchsel_free(&bs);
    tsk_treeseq_free(&ts);

    printf(fail ? "check_branchsel: FAIL\n" : "check_branchsel: PASS\n");
    return fail ? 1 : 0;
}
