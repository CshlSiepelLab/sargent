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

#include "sargent.h"
#include "branchsel.h"

/*
 * Native out-of-sample prediction (scope 3b of native_prediction_proposal.md).
 *
 * Reproduces analysis/predict.py's full five-arm ladder (sargent_s, snp_lmm,
 * sargent_a, sargent, sargent) inside the sargent binary and writes
 * <out_json_path>.  Called after a normal train fit (Stage 1 β̂ + Stage 2 ŵ²).
 *
 * Preconditions:
 *   - sg fitted on the TRAIN ARG/genotypes; sg->gr is a BED reader with
 *     precompute_stats done; sg->ws2 non-NULL (Stage 2 ran).  sg->beta may be
 *     NULL (e.g. --arg-only), treated as the zero vector.
 *   - The full ARG at full_trees_path has individuals [0,n_train) = train and
 *     [n_train,n_total) = test, matching simulate.py's layout.
 *   - Train and test BED share the same .bim (same variants, same order).
 *
 * Args:
 *   sg               fitted model (borrowed).
 *   y_train          UNCENTERED train phenotype [sg->n] (sg->y is centered).
 *   case_prefix      value written to the JSON "case" field.
 *   full_trees_path  <case>.full.trees.
 *   test_geno_prefix test BED prefix (<case>.test); .true_components.tsv is
 *                    read from <test_geno_prefix>.true_components.tsv if present.
 *   test_pheno_path  <case>.test.phen.
 *   out_json_path    output JSON (<outroot>.predict.json).
 *   bsel             branch-selection options, or NULL to skip the
 *                    sargent_w / sargent_wf arms entirely.
 *   silent           suppress progress on stderr.
 *
 * Returns 0 on success, -1 on error.
 */
/* Stage-1 CROSS-FITTING for the Stage-2 residual.  With k > 1, the training
 * PRS that forms rtrain = y - prs_train is assembled out of fold: Stage 1 is
 * refitted k times and each individual scored by a model that did not see them.
 * Without it the residual is in-sample, and a well-fitting Stage 1 drives the
 * ARG variance component to zero.  prs_TEST is unaffected, so every arm that
 * does not use rtrain is bit-identical either way.  k <= 1 disables. */
void predict_set_stage1_cv(int k);

/* --stack: emit sargent_stack, an out-of-fold-weighted blend of the Stage-1 PRS
 * and the uniform ARG BLUP.  Addresses the case where sequential
 * residualisation yields nothing because the two predictors capture the same
 * signal. */
void predict_set_stack(int on);

/* Tell predict() that beta came from a JOINT fit, so the sargent BLUP uses the
 * jointly estimated variance components instead of re-fitting REML on the
 * residual (which would discard the joint fit and reproduce the collapse). */
void predict_set_joint(int on);

int sargent_predict(sargent_t *sg,
                    const double *y_train,
                    const char *case_prefix,
                    const char *full_trees_path,
                    const char *test_geno_prefix,
                    const char *test_pheno_path,
                    const char *out_json_path,
                    const branchsel_opts_t *bsel,
                    int silent);


/* ---- Named model selection (--model) --------------------------------------
 * Historically every invocation computed every arm, so a user who wanted one
 * flavour paid for all of them: sargent_a came out identical under --snp-only,
 * --uniform-arg and --joint, and on a cold .gall cache 45-65% of the runtime
 * was work the requested model never used.  These bits let the driver skip the
 * prerequisites no selected arm needs.  SG_MODEL_ALL reproduces the old
 * behaviour exactly and stays the default. */
#define SG_ARM_SNP_LMM     (1u << 0)   /* snp_lmm     typed-SNP GRM BLUP      */
#define SG_ARM_SARGENT_S   (1u << 1)   /* sargent_s      fixed effects only      */
#define SG_ARM_SARGENT_A  (1u << 2)   /* sargent_a     uniform ARG BLUP        */
#define SG_ARM_SARGENT   (1u << 3)   /* sargent   fixed + ARG, jointly    */
/* RETIRED 2026-08-31.  Bit 4 stays reserved so no later arm silently inherits
 * a mask that older results were written under.  SARGENT-R was the Stage-2
 * learned-ws^2 predictor; it was the single best arm in 13.2% of 6452 cells,
 * its median gain over SARGENT-U was exactly 0.0000, and in the cross-ancestry
 * tail where it did win, SARGENT-UR won by more.  Localization moved to the
 * per-branch score test (Fig 5), which never used it.  Adam, 2026-08-31:
 * "the ws^2 values have not been very useful and we can retire sargent-r." */
#define SG_ARM_SARGENT_R_RETIRED (1u << 4)
#define SG_ARM_SARGENT_W   (1u << 5)   /* sargent_w   branch-selected ARG     */
#define SG_ARM_SARGENT_WF  (1u << 6)   /* sargent_wf  branch-selected + fixed */
#define SG_ARM_STACK       (1u << 7)   /* sargent_stack                       */
#define SG_MODEL_ALL       (0xFFu & ~SG_ARM_SARGENT_R_RETIRED)

/* Which arms to compute and emit.  Default SG_MODEL_ALL. */
void predict_set_models(unsigned mask);

/* Per-phase wall clock (seconds).  Set by the driver for work that happens
 * before predict_run, so one run can report a per-flavour cost instead of a
 * single process time covering every arm. */
extern double sg_t_stage1, sg_t_stage1_f, sg_t_gall, sg_t_joint, sg_t_stage2;
double arm_now(void);
