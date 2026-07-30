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
#include <string.h>
#include <getopt.h>
#include <sys/stat.h>
#include <tskit.h>
#include "geno.h"
#include "sargent.h"
#include "gall.h"
#include "variational.h"
#include "sargent.help.h"

/* Read n phenotype values from f.
 *
 * Each non-comment line is split on whitespace.  When pheno_col > 0 the
 * value at that 1-indexed field is used.  When pheno_col == 0 (legacy /
 * default), the *last* whitespace-separated field is used (single-column
 * files behave as before).
 */
static double *read_phenotypes(FILE *f, int n, int pheno_col)
{
    double *y = malloc((size_t)n * sizeof(double));
    if (!y) return NULL;

    char line[4096];
    int i = 0;
    while (i < n && fgets(line, sizeof(line), f)) {
        if (line[0] == '#') continue;
        const char *pick = NULL;
        char *tok = strtok(line, " \t\n\r");
        if (pheno_col > 0) {
            int idx = 1;
            while (tok && idx < pheno_col) { tok = strtok(NULL, " \t\n\r"); idx++; }
            if (!tok) {
                fprintf(stderr, "ERROR: phenotype line %d has fewer than %d fields\n",
                        i + 1, pheno_col);
                free(y); return NULL;
            }
            pick = tok;
        } else {
            while (tok) { pick = tok; tok = strtok(NULL, " \t\n\r"); }
            if (!pick) continue;
        }
        char *end;
        y[i] = strtod(pick, &end);
        if (end == pick) {
            fprintf(stderr, "ERROR: bad phenotype value '%s' on line %d\n",
                    pick, i + 1);
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

static void write_beta_table(sargent_t *sg, const char *outroot, int silent)
{
    int p = sg->p;
    int K = sg->K;

    char fname[4096];
    snprintf(fname, sizeof(fname), "%s.beta.tsv", outroot);
    FILE *f = fopen(fname, "w");
    if (!f) {
        fprintf(stderr, "ERROR: cannot open %s for writing\n", fname);
        return;
    }

    double *gamma_all = var_fixed_posteriors(sg);
    if (!gamma_all) {
        fprintf(stderr, "WARNING: could not compute posteriors; skipping %s\n", fname);
        fclose(f);
        return;
    }

    for (int k = 0; k < K; k++)
        fprintf(f, "gamma_%d\t", k + 1);
    fprintf(f, "E_beta\n");

    for (int j = 0; j < p; j++) {
        for (int k = 0; k < K; k++)
            fprintf(f, "%.6g\t", gamma_all[j * K + k]);
        fprintf(f, "%.6g\n", sg->beta[j]);
    }

    free(gamma_all);
    fclose(f);
    if (!silent) fprintf(stderr, "Wrote %s\n", fname);
}

static void write_w2_table(sargent_t *sg, const char *outroot, int silent)
{
    int S = sg->S;
    int M = sg->M;

    char fname[4096];
    snprintf(fname, sizeof(fname), "%s.w2.tsv", outroot);
    FILE *f = fopen(fname, "w");
    if (!f) {
        fprintf(stderr, "ERROR: cannot open %s for writing\n", fname);
        return;
    }

    for (int m = 0; m < M; m++)
        fprintf(f, "rho_%d\t", m + 1);
    fprintf(f, "E_ws2\n");

    for (int s = 0; s < S; s++) {
        for (int m = 0; m < M; m++)
            fprintf(f, "%.6g\t", sg->rho[s * M + m]);
        fprintf(f, "%.6g\n", sg->ws2[s]);
    }

    fclose(f);
    if (!silent) fprintf(stderr, "Wrote %s\n", fname);
}

/* write_summary_json: write a coherent record of the fit to <outroot>.summary.json.
 *
 * All h² values come from sargent's own internally consistent residual-based
 * calculation (h2_snp = 1 - var(r)/var(y), h2_arg = var(r)/var(y) - σ²/var(y))
 * from the final iteration's state.  PR averaging is no longer applied in
 * var_joint, so (β, ws², σ², h²) are from a single coherent iteration.
 *
 * INTERPRETATION:
 *   h2_snp = fraction of variance explained by the BayesR fit on typed SNPs.
 *   h2_arg = ADDITIONAL fraction captured by the ARG-based random component
 *            after the SNP fit ("missing heritability recovered").
 *   h2     = h2_snp + h2_arg.
 *
 * A classical fixed-vs-random variance partition is not identifiable in this
 * model — X and the ARG share structure by construction.  The
 * sequential-prioritisation interpretation above IS well-defined: it is the
 * answer the algorithm gives, not an estimate of an inaccessible truth.
 */
static void write_summary_json(sargent_t *sg, const char *outroot,
                               const char *mode, int em_iters, int silent)
{
    char fname[4096];
    snprintf(fname, sizeof(fname), "%s.summary.json", outroot);
    FILE *f = fopen(fname, "w");
    if (!f) {
        fprintf(stderr, "ERROR: cannot open %s for writing\n", fname);
        return;
    }

    fprintf(f, "{\n");
    fprintf(f, "  \"mode\":            \"%s\",\n", mode);
    fprintf(f, "  \"em_iters\":        %d,\n",     em_iters);
    fprintf(f, "  \"uniform_arg\":  %s,\n",     sg->uniform_arg ? "true" : "false");
    fprintf(f, "  \"n_indiv\":         %d,\n",     sg->n);
    fprintf(f, "  \"n_variants\":      %d,\n",     sg->p);
    fprintf(f, "  \"n_trees\":         %d,\n",     sg->S);
    fprintf(f, "  \"K_mix\":           %d,\n",     sg->K);
    fprintf(f, "  \"M_mix\":           %d,\n",     sg->M);
    fprintf(f, "  \"nprobe\":          %d,\n",     sg->nprobe);
    fprintf(f, "  \"phi_alpha\":       %.6g,\n",   sg->phi_alpha);
    fprintf(f, "  \"phi_beta\":        %.6g,\n",   sg->phi_beta);
    fprintf(f, "  \"h2\":              %.6g,\n",   sg->h2);
    fprintf(f, "  \"h2_snp\":          %.6g,\n",   sg->h2_snp);
    fprintf(f, "  \"h2_arg\":          %.6g,\n",   sg->h2_arg);
    fprintf(f, "  \"sigma2\":          %.6g\n",    sg->sigma2);
    fprintf(f, "}\n");
    fclose(f);

    if (!silent) fprintf(stderr, "Wrote %s\n", fname);
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
    int c;             /* getopt_long return; widened from signed char to int
                          to support long-only option codes >= 0x100. */
    int opt_idx;
    int silent      = 0;
    int ts_geno     = 0;
    int snp_only    = 0;
    int arg_only    = 0;
    int nprobe      = 0;   /* 0 = use sargent_t default (10) */
    int uniform_arg = 0;
    int em_iters    = 0;   /* default 0 (sequential Stage 1 then Stage 2);
                              joint EM is opt-in via --em N (N > 0) */
    int sequential_flag = 0; /* --sequential: explicit-no-op kept for forward-compat */
    int em_flag_set = 0;     /* tracks whether --em was explicitly passed */
    int pheno_col   = 0;     /* 0 = last whitespace field (legacy); >0 = 1-indexed column */
    int precompute  = 0;     /* --precompute: force rebuild of .gall eigen cache */
    double pi0      = -1.0; /* < 0 means unset */
    /* Block 2.6 — sparse fused-lasso prior on log ws². Both default to 0 (off).
     * --ws2-l1 LAMBDA   — sparsity-on-levels strength.
     * --ws2-l2 LAMBDA   — sparsity-on-differences strength. */
    double ws2_l1   = 0.0;
    double ws2_l2   = 0.0;
    const char *outroot = "sargent";
    const char *mode_str = "predict";   /* predict | localize (stub) */

    struct option long_opts[] = {
        {"silent",         0, 0, 'z'},
        {"ts-geno",        0, 0, 'G'},
        {"snp-only",       0, 0, 'S'},
        {"arg-only",       0, 0, 'A'},
        {"pi0",            1, 0, 'P'},
        {"outroot",        1, 0, 'o'},
        {"nprobe",         1, 0, 'N'},
        {"uniform-arg",    0, 0, 'U'},
        {"em",             1, 0, 'E'},
        {"sequential",     0, 0, 'Q'},
        {"pheno-col",      1, 0, 'C'},
        {"mode",           1, 0, 'M'},
        {"ws2-l1",         1, 0, 0x100},
        {"ws2-l2",         1, 0, 0x101},
        {"precompute",     0, 0, 'X'},
        {"version",        0, 0, 'v'},
        {"help",           0, 0, 'h'},
        {0, 0, 0, 0}
    };

    while ((c = getopt_long(argc, argv, "zGSAUQXE:C:M:o:N:P:vh", long_opts, &opt_idx)) != -1) {
        switch (c) {
        case 'z': silent      = 1;                break;
        case 'G': ts_geno     = 1;                break;
        case 'S': snp_only    = 1;                break;
        case 'A': arg_only    = 1;                break;
        case 'U': uniform_arg = 1;                break;
        case 'E': em_iters    = atoi(optarg);
                  em_flag_set = 1;                break;
        case 'Q': sequential_flag = 1;            break;
        case 'C': pheno_col   = atoi(optarg);     break;
        case 'M': mode_str    = optarg;           break;
        case 'P': pi0         = atof(optarg);     break;
        case 'o': outroot     = optarg;           break;
        case 'N': nprobe      = atoi(optarg);     break;
        case 0x100: ws2_l1    = atof(optarg);     break;
        case 0x101: ws2_l2    = atof(optarg);     break;
        case 'X': precompute  = 1;                break;
        case 'v': printf("sargent version %s\n", SARGENT_VERSION); exit(0);
        case 'h': printf("%s", HELP); exit(0);
        case '?':
            fprintf(stderr, "Bad argument.  Try 'sargent -h'.\n");
            exit(EXIT_FAILURE);
        }
    }

    /* Resolve mode */
    int mode_predict = (strcmp(mode_str, "predict") == 0);
    int mode_localize= (strcmp(mode_str, "localize") == 0);
    if (!mode_predict && !mode_localize) {
        fprintf(stderr, "ERROR: unrecognized --mode '%s'.  Choose from "
                "predict, localize.\n", mode_str);
        exit(EXIT_FAILURE);
    }
    if (mode_localize) {
        fprintf(stderr, "ERROR: --mode localize not yet implemented "
                "(see ROADMAP.md, Phase 7).\n");
        exit(EXIT_FAILURE);
    }

    if (snp_only && arg_only) {
        fprintf(stderr, "ERROR: --snp-only and --arg-only are mutually exclusive.\n");
        exit(EXIT_FAILURE);
    }

    /* Resolve em / sequential flags.
     *
     * Default: em_iters = 0 (sequential Stage 1 then Stage 2).
     * --em N             => opt into joint EM for N outer iterations
     *                      (joint EM is currently unstable on cases with
     *                      large numbers of trees; see ROADMAP §4.5).
     * --sequential       => explicit no-op (kept for forward-compat with
     *                      scripts written when joint EM was the default).
     * --snp-only / --arg-only force sequential and conflict with --em N>0.
     * --sequential together with an explicit --em N>0 is an error.
     */
    if (sequential_flag && em_flag_set && em_iters > 0) {
        fprintf(stderr, "ERROR: --sequential conflicts with --em N (N>0).\n");
        exit(EXIT_FAILURE);
    }
    if (sequential_flag) em_iters = 0;
    if (snp_only || arg_only) {
        if (em_flag_set && em_iters > 0) {
            fprintf(stderr, "ERROR: --em cannot be combined with --snp-only or --arg-only.\n");
            exit(EXIT_FAILURE);
        }
        em_iters = 0;
    }

    int n_pos   = argc - optind;
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

    FILE *pheno_fh = fopen(pheno_file, "r");
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

    if (!silent) fprintf(stderr, "Reading %s...\n", trees_file);
    tsk_treeseq_t ts;
    CHECK(tsk_treeseq_load(&ts, trees_file, 0));

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

    if (!silent) fprintf(stderr, "Precomputing variant statistics...\n");
    if (geno_precompute_stats(gr) < 0) {
        fprintf(stderr, "ERROR: precompute_stats failed.\n");
        geno_close(gr);
        tsk_treeseq_free(&ts);
        exit(EXIT_FAILURE);
    }

    if (!silent)
        fprintf(stderr, "Polymorphic variants: %d\n", gr->n_variants);

    double *y = read_phenotypes(pheno_fh, gr->n_samples, pheno_col);
    fclose(pheno_fh);
    if (!y) {
        geno_close(gr);
        tsk_treeseq_free(&ts);
        exit(EXIT_FAILURE);
    }

    sargent_t *sg = sargent_new(gr, &ts, y);
    free(y);
    if (!sg) {
        fprintf(stderr, "ERROR: failed to build sargent object.\n");
        geno_close(gr);
        tsk_treeseq_free(&ts);
        exit(EXIT_FAILURE);
    }
    if (nprobe > 0)  sg->nprobe = nprobe;
    if (uniform_arg) sg->uniform_arg = 1;
    if (em_iters > 0) sg->em_iters = em_iters;
    if (ws2_l1 > 0.0) sg->fused_lasso_lambda1 = ws2_l1;
    if (ws2_l2 > 0.0) sg->fused_lasso_lambda2 = ws2_l2;
    if (pi0 > 0.0 && pi0 < 1.0) {
        /* Beta(1, pi0/(1-pi0)) hyperprior: prior mean of phi = 1-pi0.
         * E.g. pi0=0.9 → Beta(1,9), prior mean phi=0.1 (10% active). */
        sg->phi_alpha = 1.0;
        sg->phi_beta  = pi0 / (1.0 - pi0);
    }
    (void)mode_predict;  /* currently the only operational mode */

    /* ---- G_all precompute / cache load ------------------------------------
     *
     * G_all = Σ_s G_s enables two optimisations in Stage 2:
     *   1. EMMA-based Newton-REML for --uniform-arg: O(n²) vs. O(100 × 700 × S)
     *      for CG-based var_random (>60× speedup at n=1000, S=3000).
     *   2. G_all preconditioner for per-tree CG: κ(P⁻¹V)≈O(1) → ~5 iterations
     *      vs. ~64 without preconditioner.
     *
     * Cache file: <trees_stem>.gall beside the .trees file.  Invalidated when
     * the .trees file size changes.  Auto-detected on every run; rebuilt when
     * --precompute is passed or on first use with --uniform-arg.
     * -------------------------------------------------------------------- */
    gall_t *gall = NULL;
    {
        /* Derive .gall path from .trees path by stripping extension. */
        char gall_path[4096];
        snprintf(gall_path, sizeof(gall_path), "%s", trees_file);
        char *ext = strrchr(gall_path, '.');
        if (ext && strcmp(ext, ".trees") == 0) *ext = '\0';
        strncat(gall_path, ".gall", sizeof(gall_path) - strlen(gall_path) - 1);

        /* Get current .trees file size for cache validation. */
        struct stat ts_stat;
        long trees_sz = (stat(trees_file, &ts_stat) == 0) ? (long)ts_stat.st_size : -1L;

        if (precompute) {
            if (!silent)
                fprintf(stderr, "G_all: --precompute: rebuilding %s\n", gall_path);
        } else {
            /* Try loading cached .gall */
            gall = gall_load(gall_path, sg->n, trees_sz);
            if (gall && !silent)
                fprintf(stderr, "G_all: loaded %s\n", gall_path);
        }

        /* Build (and save) if not loaded and we need it:
         *   - always for --uniform-arg (EMMA path requires gall)
         *   - always for --precompute
         *   - for per-tree mode: gall supplies both the CG preconditioner and
         *     the EMMA warm-start (reduces Stage 2 from ~750 s to ~50 s). */
        if (!gall && (precompute || uniform_arg || !snp_only)) {
            /* Span-weight G_all only for the uniform arms (their kernel is the
             * span-weighted eGRM). For the learned arms G_all is used solely as the
             * CG preconditioner, which must match the un-span learned operator. */
            gall = gall_build(&ts, sg->hap_to_indiv, sg->n, uniform_arg);
            if (gall) {
                if (trees_sz >= 0) {
                    if (gall_save(gall, gall_path, trees_sz) == 0 && !silent)
                        fprintf(stderr, "G_all: saved %s\n", gall_path);
                }
            } else {
                fprintf(stderr, "WARNING: G_all build failed; "
                        "falling back to CG-based Stage 2\n");
            }
        }

        sg->gall = gall;
    }

    if (em_iters > 0) {
        /* ---- Warm-started joint EM (Phase 4.5) -----------------------------
         *
         * Cold-start var_joint was empirically unstable at scale: on the
         * Nurdan smoke case (n=200, p=2488, S=2739) the β-sweep with stale
         * r̃ caused β to explode in iter 1, which drove σ² oscillation
         * (collapse with the original A_0/B_0 update; bounded oscillation
         * after Phase 4.5 damping + ratio-clipping).  Damping alone did
         * not yield convergence.
         *
         * Phase 4.5 ships warm-start: run sequential (var_fixed then
         * var_random) first, then use that β/ws²/σ² as the initial point
         * for var_joint.  Treats joint EM as a polish step on a known-good
         * fit rather than as a cold-start optimiser.
         * ------------------------------------------------------------------ */
        if (!silent)
            fprintf(stderr,
                    "Joint EM: warm-start (sequential + %d polish iters; "
                    "%d SNPs K=%d, %d trees%s)\n",
                    em_iters, sg->p, sg->K, sg->S,
                    sg->uniform_arg ? ", uniform-random" : "");

        if (!silent) fprintf(stderr, "  warm-start Stage 1 (BayesR fit to SNPs)...\n");
        var_fixed(sg, silent);
        if (!silent)
            fprintf(stderr, "  warm-start Stage 2 (per-tree REML EM)...\n");
        var_random(sg, silent);

        if (!silent)
            fprintf(stderr,
                    "  joint EM polish (%d outer iters from warm β/ws²/σ²)\n",
                    em_iters);
        var_joint(sg, silent);
        if (!silent)
            fprintf(stderr, "Joint EM done: h2_snp=%.4f  h2_arg=%.4f  sigma2=%.4f\n",
                    sg->h2_snp, sg->h2_arg, sg->sigma2);
        write_beta_table(sg, outroot, silent);
        write_w2_table(sg, outroot, silent);
    } else {
        /* ---- Sequential mode: Stage 1 then Stage 2 (mode: predict) ---- */
        if (!arg_only) {
            if (!silent)
                fprintf(stderr, "Stage 1: inferring fixed effects (%d SNPs, K=%d)...\n",
                        sg->p, sg->K);
            var_fixed(sg, silent);
            if (!silent)
                fprintf(stderr, "Stage 1 done: h2_snp=%.4f  sigma2=%.4f\n",
                        sg->h2_snp, sg->sigma2);
            write_beta_table(sg, outroot, silent);
        }

        if (!snp_only) {
            int use_emma       = (sg->uniform_arg && sg->gall != NULL);
            int emma_warmstart = (!sg->uniform_arg && sg->gall != NULL);
            if (!silent)
                fprintf(stderr, "Stage 2: inferring segment weights (%d trees%s%s)...\n",
                        sg->S,
                        sg->uniform_arg ? ", uniform-random" : "",
                        use_emma ? ", EMMA eigenspace" :
                        emma_warmstart ? ", EMMA warm-start" : "");
            if (use_emma) {
                var_random_emma(sg, silent);
            } else {
                if (emma_warmstart) {
                    /* Run EMMA in uniform mode to obtain a fast initialiser for
                     * ws2[s] and sigma2.  Temporarily set uniform_arg=1 so that
                     * var_random_emma takes the EMMA path; restore before the
                     * full per-tree EM call.  On EMMA boundary failure the
                     * function falls back to uniform EM (still a better warm-start
                     * than tau2_mix/2).  var_random skips re-initialisation when
                     * ws2 is already non-NULL, so the warm-start is picked up
                     * automatically. */
                    sg->uniform_arg = 1;
                    var_random_emma(sg, /*silent=*/1);
                    sg->uniform_arg = 0;
                }
                var_random(sg, silent);
            }
            if (!silent)
                fprintf(stderr, "Stage 2 done: h2_arg=%.4f  sigma2=%.4f\n",
                        sg->h2_arg, sg->sigma2);
            write_w2_table(sg, outroot, silent);
        }
    }

    /* Write coherent summary of the final fit. */
    {
        const char *m = em_iters > 0 ? "joint_em" : "sequential";
        write_summary_json(sg, outroot, m, em_iters, silent);
    }

    /* ---- Summary statistics ---- */
    if (!silent) {
        fprintf(stderr, "\nSummary:\n");
        fprintf(stderr, "  h2        : %.4f\n", sg->h2);
        fprintf(stderr, "  h2.snp    : %.4f\n", sg->h2_snp);
        fprintf(stderr, "  h2.arg    : %.4f\n", sg->h2_arg);
        fprintf(stderr, "  sigma2    : %.4f\n", sg->sigma2);
    }

    sg->gall = NULL;       /* gall is managed externally; clear before free */
    sargent_free(sg);
    gall_free(gall);
    geno_close(gr);
    tsk_treeseq_free(&ts);
    if (!silent) fprintf(stderr, "Done.\n");
    return 0;
}
