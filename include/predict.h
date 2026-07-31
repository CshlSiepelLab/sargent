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
 *   silent           suppress progress on stderr.
 *
 * Returns 0 on success, -1 on error.
 */
int sargent_predict(sargent_t *sg,
                    const double *y_train,
                    const char *case_prefix,
                    const char *full_trees_path,
                    const char *test_geno_prefix,
                    const char *test_pheno_path,
                    const char *out_json_path,
                    int silent);
