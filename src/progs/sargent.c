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
#include "predict.h"
#include "threads.h"
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
    int stage1_gibbs = 0, gibbs_iter = 0, gibbs_burn = -1;
    int stage1_cv = 0, use_stack = 0, use_joint = 0;
    /* --model: name the flavour instead of spelling out the flag combination,
     * and run ONLY that flavour's prerequisites.  NULL = no --model given, in
     * which case the legacy flags decide and the mask stays SG_MODEL_ALL so
     * existing scripts are bit-identical. */
    const char *model_name = NULL;
    int model_skip_gall = 0;
    unsigned arm_mask = SG_MODEL_ALL;
    int arg_only    = 0;
    int nprobe      = 0;   /* 0 = use sargent_t default (10) */
    int uniform_arg = 0;
    int em_iters    = 0;   /* default 0 (sequential Stage 1 then Stage 2);
                              joint EM is opt-in via --em N (N > 0) */
    int sequential_flag = 0; /* --sequential: explicit-no-op kept for forward-compat */
    int em_flag_set = 0;     /* tracks whether --em was explicitly passed */
    int pheno_col   = 0;     /* 0 = last whitespace field (legacy); >0 = 1-indexed column */
    int precompute  = 0;     /* --precompute: force rebuild of .gall eigen cache */
    int threads     = 0;     /* --threads N; 0 = resolve from env, else library default */
    double pi0      = -1.0; /* < 0 means unset */
    /* Block 2.6 — sparse fused-lasso prior on log ws². Both default to 0 (off).
     * --ws2-l1 LAMBDA   — sparsity-on-levels strength.
     * --ws2-l2 LAMBDA   — sparsity-on-differences strength. */
    double ws2_l1   = 0.0;
    double ws2_l2   = 0.0;
    const char *outroot = "sargent";
    const char *mode_str = "predict";   /* predict | localize (stub) */

    /* Native out-of-sample prediction (--predict).  Paths default from the
     * geno_prefix (the case prefix) per simulate.py's conventions. */
    int predict_oos = 0;
    const char *full_trees_arg = NULL;   /* default <geno_prefix>.full.trees */
    const char *test_geno_arg  = NULL;   /* default <geno_prefix>.test */
    const char *test_pheno_arg = NULL;   /* default <geno_prefix>.test.phen */
    const char *predict_out_arg = NULL;  /* default <outroot>.predict.json */

    /* Branch selection (--branch-select): the Stage-2 replacement.  Adds the
     * sargent_w / sargent_wf arms to --predict.  Opt-in for now; flipping it to
     * the default is deferred until Figs 2-4 have been regenerated on it. */
    int bsel_on    = 0;
    int bsel_k     = BRANCHSEL_DEFAULT_K;
    int bsel_trees = BRANCHSEL_DEFAULT_TREES;

    struct option long_opts[] = {
        {"silent",         0, 0, 'z'},
        {"ts-geno",        0, 0, 'G'},
        {"snp-only",       0, 0, 'S'},
        {"arg-only",       0, 0, 'A'},
        {"pi0",            1, 0, 'P'},
        {"outroot",        1, 0, 'o'},
        {"nprobe",         1, 0, 'N'},
        {"threads",        1, 0, 'T'},
        {"uniform-arg",    0, 0, 'U'},
        {"em",             1, 0, 'E'},
        {"sequential",     0, 0, 'Q'},
        {"pheno-col",      1, 0, 'C'},
        {"mode",           1, 0, 'M'},
        {"ws2-l1",         1, 0, 0x100},
        {"ws2-l2",         1, 0, 0x101},
        {"predict",        0, 0, 0x102},
        {"full-trees",     1, 0, 0x103},
        {"test-geno",      1, 0, 0x104},
        {"test-pheno",     1, 0, 0x105},
        {"predict-out",    1, 0, 0x106},
        {"model",          1, 0, 0x130},
        {"stage1-gibbs",   0, 0, 0x107},
        {"stage1-cv",      1, 0, 0x120},
        {"stack",          0, 0, 0x121},
        {"joint",          0, 0, 0x122},
        {"gibbs-iter",     1, 0, 0x108},
        {"gibbs-burn",     1, 0, 0x109},
        {"branch-select",  0, 0, 0x10a},
        {"select-k",       1, 0, 0x10b},
        {"select-trees",   1, 0, 0x10c},
        {"precompute",     0, 0, 'X'},
        {"version",        0, 0, 'v'},
        {"help",           0, 0, 'h'},
        {0, 0, 0, 0}
    };

    while ((c = getopt_long(argc, argv, "zGSAUQXE:C:M:o:N:P:T:vh", long_opts, &opt_idx)) != -1) {
        switch (c) {
        case 'z': silent      = 1;                break;
        case 'G': ts_geno     = 1;                break;
        case 'S': snp_only    = 1;                break;
        case 0x130: model_name = optarg;          break;
        case 0x107: stage1_gibbs = 1;             break;
        case 0x120: stage1_cv = atoi(optarg);    break;
        case 0x121: use_stack = 1;               break;
        case 0x122: use_joint = 1;               break;
        case 0x108: gibbs_iter   = atoi(optarg);  break;
        case 0x109: gibbs_burn   = atoi(optarg);  break;
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
        case 'T': threads     = atoi(optarg);     break;
        case 0x100: ws2_l1    = atof(optarg);     break;
        case 0x101: ws2_l2    = atof(optarg);     break;
        case 0x102: predict_oos    = 1;           break;
        case 0x103: full_trees_arg = optarg;      break;
        case 0x104: test_geno_arg  = optarg;      break;
        case 0x105: test_pheno_arg = optarg;      break;
        case 0x106: predict_out_arg= optarg;      break;
        case 0x10a: bsel_on    = 1;               break;
        case 0x10b: bsel_k     = atoi(optarg); bsel_on = 1; break;
        case 0x10c: bsel_trees = atoi(optarg); bsel_on = 1; break;
        case 'X': precompute  = 1;                break;
        case 'v': printf("sargent version %s\n", SARGENT_VERSION); exit(0);
        case 'h': printf("%s", HELP); exit(0);
        case '?':
            fprintf(stderr, "Bad argument.  Try 'sargent -h'.\n");
            exit(EXIT_FAILURE);
        }
    }

    /* Establish one thread budget across OpenMP and BLAS before any linear
     * algebra runs.  See include/threads.h for why a single env var was not
     * enough. */
    sargent_init_threads(threads, silent);

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

    /* ---- --model: name -> arms -> prerequisites ---------------------------
     * Each model implies which arms are emitted AND which expensive stages have
     * to run.  Skipping the rest is the point: on a cold .gall cache the
     * Stage-1 Gibbs is ~186 s and building G_all ~123 s, so a user who wants
     * only SARGENT-UR was paying ~186 s of Stage 1 they never used, and one who
     * wants only SARGENT-F ~123 s of genealogy work they never used. */
    if (model_name) {
        struct { const char *name; unsigned arms; } models[] = {
            { "snp-lmm",    SG_ARM_SNP_LMM    },
            { "sargent-f",  SG_ARM_SARGENT_F  },
            { "sargent-ur", SG_ARM_SARGENT_UR },
            { "sargent-u",  SG_ARM_SARGENT_U  },
            { "sargent-r",  SG_ARM_SARGENT_R  },
            { "sargent-w",  SG_ARM_SARGENT_W  },
            { "sargent-wf", SG_ARM_SARGENT_WF },
            /* "all" is the FIGURE set, and deliberately excludes SARGENT-R.
             * SARGENT-R is a different Stage 2 -- the per-tree learned ws^2 fit
             * -- and it is ~19x the cost of every other arm combined (3992 s
             * against 215 s on the benchmark cell).  Bundling it into "all"
             * made --model all 9x slower than the invocation it was meant to
             * replace.  Excluding it also makes "all" an exact drop-in for what
             * every figure driver already runs (--stage1-gibbs --joint
             * --branch-select --uniform-arg); SARGENT-R stays available as its
             * own model, which is how the repo already harvests it (into
             * est_r.predict.json, via a separate --arg-only pass). */
            { "all",        SG_MODEL_ALL & ~SG_ARM_SARGENT_R },
        };
        size_t nm = sizeof(models) / sizeof(models[0]), mi;
        for (mi = 0; mi < nm; mi++)
            if (strcmp(model_name, models[mi].name) == 0) break;
        if (mi == nm) {
            fprintf(stderr, "ERROR: unknown --model '%s'.  Choose one of:", model_name);
            for (size_t k = 0; k < nm; k++) fprintf(stderr, " %s", models[k].name);
            fprintf(stderr, "\n");
            return 1;
        }
        arm_mask = models[mi].arms;

        /* Derive the PREREQUISITES from the arms, not from the legacy flags.
         * Overloading --snp-only / --arg-only here was wrong: they are a blunt
         * pair (one skips Stage 2, the other skips Stage 1) and setting both, or
         * setting --arg-only alongside --stage1-gibbs, trips the parser's own
         * consistency checks.  These booleans say exactly what each model needs. */
        /* SARGENT-R is the ARG-ONLY learned-ws^2 model: the repo harvests it with
         * `--arg-only --predict` and reads the `sargent` key, so it needs Stage 2
         * but NOT Stage 1.  Listing it under need_stage1 was wrong. */
        int need_stage1 = (arm_mask & (SG_ARM_SARGENT_F | SG_ARM_SARGENT_U |
                                       SG_ARM_SARGENT_WF)) != 0;
        int need_gall   = (arm_mask & (SG_ARM_SARGENT_UR | SG_ARM_SARGENT_U |
                                       SG_ARM_SARGENT_R  | SG_ARM_SARGENT_W |
                                       SG_ARM_SARGENT_WF)) != 0;
        int need_stage2 = (arm_mask & SG_ARM_SARGENT_R) != 0;

        if (arm_mask & (SG_ARM_SARGENT_U | SG_ARM_SARGENT_WF)) use_joint = 1;
        if (arm_mask & (SG_ARM_SARGENT_W | SG_ARM_SARGENT_WF)) bsel_on   = 1;
        if (!need_stage2) uniform_arg = 1;   /* only sargent-r fits the per-tree ws^2 */
        if (!need_stage1) arg_only    = 1;   /* --arg-only IS "skip Stage 1"          */
        else              stage1_gibbs = 1;  /* variational Stage 1 is not at parity  */
        model_skip_gall = !need_gall;        /* checked at the .gall build site       */
    }
    predict_set_models(arm_mask);

    var_fixed_set_gibbs(stage1_gibbs, gibbs_iter, gibbs_burn);
    predict_set_stage1_cv(stage1_cv);
    predict_set_stack(use_stack);
    predict_set_joint(use_joint);
    if (stage1_gibbs && arg_only) {
        fprintf(stderr, "ERROR: --stage1-gibbs has no effect with --arg-only "
                        "(Stage 1 is skipped).\n");
        return 1;
    }

    if (snp_only && arg_only) {
        fprintf(stderr, "ERROR: --snp-only and --arg-only are mutually exclusive.\n");
        exit(EXIT_FAILURE);
    }

    if (bsel_on) {
        if (!predict_oos) {
            fprintf(stderr, "ERROR: --branch-select only has an effect with "
                    "--predict (it produces prediction arms, not a fit).\n");
            exit(EXIT_FAILURE);
        }
        if (bsel_k < 0) {
            fprintf(stderr, "ERROR: --select-k must be >= 0.\n");
            exit(EXIT_FAILURE);
        }
        if (bsel_trees < 0) {
            fprintf(stderr, "ERROR: --select-trees must be >= 0 (0 = all trees).\n");
            exit(EXIT_FAILURE);
        }
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

    if (predict_oos) {
        if (ts_geno) {
            fprintf(stderr, "ERROR: --predict requires BED genotypes "
                    "(not --ts-geno); test genotypes come from a test BED.\n");
            exit(EXIT_FAILURE);
        }
        /* --snp-only is ALLOWED.  Four of the five arms (bayesR, snp_lmm,
         * arg_lmm, sargent_u) never touch the learned per-tree ws2 -- the
         * uniform kinship accumulator in predict.c is independent of it.  Only
         * the full `sargent` arm needs Stage 2, and predict_oos() omits that arm
         * from the JSON when ws2 is absent.  Stage 2 is ~88% of a cell, so this
         * is the cheap path for sweeps that do not need the learned-weight arm. */
        if (snp_only) {
            fprintf(stderr, "note: --predict with --snp-only reports four arms "
                    "(bayesR, snp_lmm, arg_lmm, sargent_u); `sargent` needs the "
                    "Stage-2 ws2 fit and is omitted.\n");
        }
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
    /* Keep the uncentered `y` alive: --predict needs the original train
     * phenotype (sg->y is centered).  Freed at the end of main(). */
    if (!sg) {
        free(y);
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
    /* A model with no ARG arm (snp-lmm, sargent-f) must not pay for G_all:
     * building and eigendecomposing it is ~123 s on a 500-sample cell. */
    int skip_gall = model_skip_gall;
    if (!skip_gall) {
        /* Derive the .gall path from the .trees path by stripping the extension.
         *
         * TWO KERNELS, TWO FILES.  gall_build produces the span-weighted eGRM
         * (Σ_s span_s·G_s) under --uniform-arg and the plain Σ_s G_s otherwise;
         * these are different matrices and, for --uniform-arg, G_all is the
         * MODEL rather than a preconditioner.  Before v4 both shared one path,
         * so the driver order (`--predict` → `--arg-only` → `--uniform-arg`,
         * all in one cell dir) had the uniform pass silently adopt the learned
         * pass's un-span cache.  Separate names keep BOTH caches valid, so
         * neither pass has to rebuild the other's; the use_span field in the
         * header (checked by gall_load) is the belt-and-braces guard for any
         * pre-existing or hand-copied file. */
        char gall_path[4096];
        snprintf(gall_path, sizeof(gall_path), "%s", trees_file);
        char *ext = strrchr(gall_path, '.');
        if (ext && strcmp(ext, ".trees") == 0) *ext = '\0';
        /* ONE kernel, one file: span weighting is unconditional since the
         * per-segment weights were retired, so the .span.gall / .gall split is
         * gone.  The v5 header rejects any file written under the old scheme. */
        strncat(gall_path, ".gall",
                sizeof(gall_path) - strlen(gall_path) - 1);

        /* Get current .trees file size for cache validation. */
        struct stat ts_stat;
        long trees_sz = (stat(trees_file, &ts_stat) == 0) ? (long)ts_stat.st_size : -1L;

        if (precompute) {
            if (!silent)
                fprintf(stderr, "G_all: --precompute: rebuilding %s\n", gall_path);
        } else {
            /* Try loading cached .gall (must hold the kernel this run needs) */
            { double _t = arm_now(); gall = gall_load(gall_path, sg->n, trees_sz, /*use_span=*/1); sg_t_gall += arm_now() - _t; }
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
            { double _t = arm_now(); gall = gall_build(&ts, sg->hap_to_indiv, sg->n, /*use_span=*/1); sg_t_gall += arm_now() - _t; }
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
        { double _t = arm_now(); var_fixed(sg, silent); sg_t_stage1 += arm_now() - _t; }
        if (!silent)
            fprintf(stderr, "  warm-start Stage 2 (per-tree REML EM)...\n");
        { double _t = arm_now(); var_random(sg, silent); sg_t_stage2 += arm_now() - _t; }

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
            /* --joint replaces the SEQUENTIAL Stage 1 + Stage 2 with a single
             * chain that fits beta and the ARG variance components together.
             * Sequential gives the ARG nothing -- sargent_u equalled bayesR to
             * six decimals at every phi level -- because beta gets first claim
             * on everything the two terms share.  See var_joint_gibbs. */
            if (use_joint) {
                if (!sg->gall) {
                    fprintf(stderr, "ERROR: --joint needs the .gall eigenbasis; "
                                    "it is unavailable here (per-tree kernel or "
                                    "cache miss)\n");
                    return 1;
                }
                /* SARGENT-F is the fixed-effects-only model and has no ARG
                 * term, so there is nothing in it to estimate jointly: its beta
                 * must come from fitting y = X.beta + eps on its own.  Do that
                 * FIRST and stash the result, then let the joint fit overwrite
                 * sg->beta.  Costs one extra Stage 1, which is the honest price
                 * of reporting two different models from one invocation. */
                if (!silent)
                    fprintf(stderr, "  standalone Stage 1 for SARGENT-F "
                                    "(fixed effects only)...\n");
                { double _t = arm_now(); { double _t = arm_now(); var_fixed(sg, silent); sg_t_stage1 += arm_now() - _t; }
                  sg_t_stage1_f = arm_now() - _t; }
                sg->beta_f = malloc((size_t)sg->p * sizeof(double));
                if (!sg->beta_f) { fprintf(stderr, "ERROR: out of memory\n"); return 1; }
                memcpy(sg->beta_f, sg->beta, (size_t)sg->p * sizeof(double));

                if (!silent)
                    fprintf(stderr, "Stage 1+2 JOINT: beta and ARG variance "
                                    "components in one chain (%d SNPs, K=%d)...\n",
                            sg->p, sg->K);
                double j_su2 = 0.0, j_sig2 = 0.0;
                double _tj = arm_now();
                if (var_joint_gibbs(sg, sg->gall, silent, gibbs_iter ? gibbs_iter : 3000,
                                    gibbs_burn >= 0 ? gibbs_burn : 1000,
                                    0, &j_su2, &j_sig2) != 0) {
                    fprintf(stderr, "ERROR: joint estimation failed\n");
                    return 1;
                }
                sg->joint_su2 = j_su2;
                sg->sigma2    = j_sig2;
                if (!silent)
                    fprintf(stderr, "Joint done: su2=%.4f sigma2=%.4f\n", j_su2, j_sig2);
                write_beta_table(sg, outroot, silent);
            } else {
            if (!silent)
                fprintf(stderr, "Stage 1: inferring fixed effects (%d SNPs, K=%d)...\n",
                        sg->p, sg->K);
            { double _t = arm_now(); var_fixed(sg, silent); sg_t_stage1 += arm_now() - _t; }
            if (!silent)
                fprintf(stderr, "Stage 1 done: h2_snp=%.4f  sigma2=%.4f\n",
                        sg->h2_snp, sg->sigma2);
            write_beta_table(sg, outroot, silent);
            }
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
                { double _t = arm_now(); var_random_emma(sg, silent); sg_t_stage2 += arm_now() - _t; }
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
                { double _t = arm_now(); var_random(sg, silent); sg_t_stage2 += arm_now() - _t; }
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

    /* ---- Native out-of-sample prediction (--predict) --------------------- */
    int predict_rc = 0;
    if (predict_oos) {
        char full_trees[4096], test_geno[4096], test_pheno[4096], predict_out[4096];
        if (full_trees_arg) snprintf(full_trees, sizeof(full_trees), "%s", full_trees_arg);
        else                snprintf(full_trees, sizeof(full_trees), "%s.full.trees", geno_prefix);
        if (test_geno_arg)  snprintf(test_geno, sizeof(test_geno), "%s", test_geno_arg);
        else                snprintf(test_geno, sizeof(test_geno), "%s.test", geno_prefix);
        if (test_pheno_arg) snprintf(test_pheno, sizeof(test_pheno), "%s", test_pheno_arg);
        else                snprintf(test_pheno, sizeof(test_pheno), "%s.test.phen", geno_prefix);
        if (predict_out_arg) snprintf(predict_out, sizeof(predict_out), "%s", predict_out_arg);
        else                 snprintf(predict_out, sizeof(predict_out), "%s.predict.json", outroot);

        if (!silent) fprintf(stderr, "\nStage 3: out-of-sample prediction ...\n");
        branchsel_opts_t bopts = { .k = bsel_k, .max_trees = bsel_trees,
                                   .silent = silent };
        predict_rc = sargent_predict(sg, y, geno_prefix, full_trees,
                                     test_geno, test_pheno, predict_out,
                                     bsel_on ? &bopts : NULL, silent);
        if (predict_rc != 0)
            fprintf(stderr, "ERROR: prediction stage failed.\n");
    }

    sg->gall = NULL;       /* gall is managed externally; clear before free */
    sargent_free(sg);
    free(y);
    gall_free(gall);
    geno_close(gr);
    tsk_treeseq_free(&ts);
    if (!silent) fprintf(stderr, "Done.\n");
    return predict_rc == 0 ? 0 : EXIT_FAILURE;
}
