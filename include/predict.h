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
 * Reproduces analysis/predict.py's full five-arm ladder (bayesR, snp_lmm,
 * arg_lmm, sargent_u, sargent) inside the sargent binary and writes
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

/* Tell predict() that beta came from a JOINT fit, so the sargent_u BLUP uses the
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

/* Per-phase wall clock (seconds).  Set by the driver for work that happens
 * before predict_run, so one run can report a per-flavour cost instead of a
 * single process time covering every arm. */
extern double sg_t_stage1, sg_t_stage1_f, sg_t_gall, sg_t_joint, sg_t_stage2;
double arm_now(void);
