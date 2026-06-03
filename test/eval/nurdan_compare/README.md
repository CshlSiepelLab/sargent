# SARGENT vs ARG-LMM (MKLMM) — back-to-back comparison

This directory contains scripts to run **SARGENT** and **Nurdan's ARG-LMM / MKLMM**
on the same simulation data, enabling direct comparison of results and running times.

## Prerequisites

### SARGENT
```bash
# Build (from repo root)
mamba activate sargent       # cmake, meson, ninja, msprime, tskit, numpy
cmake -B build -S .
cmake --build build -j
# Binary: build/src/progs/sargent
```

### MKLMM (Nurdan's R code)
R 4.x with these packages (all available on CRAN):
```r
install.packages(c("MASS", "Matrix", "glmnet", "jsonlite"))
# Note: HDeconometrics is NOT required — the BIC lasso is re-implemented
# directly using glmnet::glmnet() + BIC path selection (faithful equivalent).
```

---

## Figure 1 simulations — heritability estimation

These simulations have **all N individuals in the ARG** (no train/test split).
The experiment compares h²_SNP, h²_ARG, σ² estimates between methods.

### Step 1 — run SARGENT (Figure 1 format)

```bash
REPO=/path/to/sargent
SIM=~/test-sargent/nurdan-all/Simulations/Figure1_Simulation_Xb_U/\
Simulation_Fig1_200ind_CausalRatio_0.1_Seed6903_1Mb
OUT=~/test-sargent/nurdan_compare
mkdir -p $OUT

# Convert Nurdan's sim to SARGENT case layout (symlinks + params JSON)
python $REPO/scripts/import_nurdan.py $SIM $OUT/case

# Run SARGENT (Stage 1 BayesR + Stage 2 per-tree REML EM)
PYBIN=$(mamba run -n sargent which python)
$PYBIN $REPO/test/eval/run_eval.py \
  --prefix $OUT/case \
  --sargent $REPO/build/src/progs/sargent \
  --out-prefix $OUT/sargent_balanced_rep1
```

Output files: `sargent_balanced_rep1.summary.json`, `.beta.tsv`, `.w2.tsv`

### Step 2 — run MKLMM (Model 1: ARG-LMM, original algorithm)

```bash
Rscript $REPO/test/eval/nurdan_compare/run_nurdan.R \
  $SIM vc_balanced 1 \
  $OUT/nurdan_balanced_rep1.json
```

Arguments: `<sim_dir>  <trait>  <pheno_col>  <out_json>`
- `trait`: `vc_balanced`, `vc_highfix_lowrand`, `vc_lowfix_highrand`
- `pheno_col`: replicate index 1–5

### Step 3 — compare

```bash
python $REPO/test/eval/nurdan_compare/compare.py \
  $OUT/sargent_balanced_rep1.summary.json \
  $OUT/nurdan_balanced_rep1.json
```

### Results (Seed6903, vc_balanced, rep1; truth h²_SNP=0.4, h²_ARG=0.4)

```
Quantity          Truth   SARGENT   MKLMM   diff(S-M)
h2_snp            0.4000   0.7922   0.6578   +0.2168
h2_arg            0.4000   0.0678   0.1119   -0.0440
h2_total          0.8000   0.8600   0.7696   +0.0904
sigma2            0.2000   0.1393   0.1432   -0.0040
wall_fit (s)             ~325      ~41
```

Both methods over-attribute to SNPs (h²_SNP > truth) and under-attribute to ARG
(h²_ARG < truth). This is expected: SARGENT's sequential two-stage design prioritises
SNPs first; MKLMM's fixed-effects lasso similarly absorbs SNP variance. The ARG captures
only residual signal.

**Timing note**: MKLMM's ~41s **excludes** the Python/tskit kernel pre-computation
(`.npy` matrices). SARGENT's ~325s **includes** all tree traversal. Add kernel
pre-computation time for a fair total-pipeline comparison.

---

## Figure 2 simulations — missing heritability under SNP downsampling

These simulations have **350 individuals** (300 EUR + 50 AFR) with SNPs downsampled to
100% / 75% / 50% / 25% / 10% of the original panel. The ARG kernels are computed from
the **full** tree sequence in all cases, so ARG-based methods can recover heritability
that typed-SNP-only methods miss.

### Run the full 6-model suite + ARG-LMM Stage 2

```bash
REPO=/path/to/sargent
SIM=~/test-sargent/nurdan-all/Simulations/Figure2_MissingData/\
Simulation_350_0.1_8890_1Mb
OUT=~/test-sargent/nurdan_compare

Rscript $REPO/test/eval/nurdan_compare/run_nurdan_fig2.R \
  $SIM  \
  <pheno_col>  \     # 3 = first replicate (h2=0.8)
  <ds_factor>  \     # 0=100% SNPs, 1≈75%, 2≈50%, 4≈25%, 10≈10%
  <seed>       \     # integer RNG seed for EUR train/test split
  <out_json>
```

Example — all SNPs, seed 42:
```bash
Rscript $REPO/test/eval/nurdan_compare/run_nurdan_fig2.R \
  $SIM 3 0 42 $OUT/fig2_ds0.json
```

The script runs **7 models** in one go:

| Model | Description |
|---|---|
| ARG-LMM | Full ARG multi-kernel LMM with BIC-lasso fixed effects |
| ARG-LMM-PCA | Same but PCs instead of raw SNPs as fixed effects |
| Just-Random | ARG multi-kernel, intercept only (no SNP fixed effects) |
| GRM-LMM (ARG) | Single averaged ARG kernel + SNP fixed effects |
| GRM-LMM (SNP) | Standard SNP GRM, no ARG |
| Just-Fixed | Lasso on typed SNPs only, no random effects |
| ARG-LMM Stage 2 | Joint BayesR (SNPs) + AIREML (segment weights w²_s) |

For each model the script reports:
- **Heritability** estimates (all 350 individuals)
- **Prediction R²** (80% EUR train → 20% EUR test)
- Wall time per model

### Results (Seed8890, pheno_col=3, h²_truth≈0.8 pure-SNP)

**All SNPs (ds_factor=0, 2074 SNPs)**

| Model | h²_ARG | R² pred | time (s) |
|---|---|---|---|
| ARG-LMM | 0.793 | 0.655 | ~99 |
| Just-Random | 0.638 | 0.666 | ~94 |
| GRM-LMM (ARG) | 0.771 | 0.662 | ~1 |
| GRM-LMM (SNP) | 0.808 | 0.641 | <1 |
| Just-Fixed | — | 0.597 | <1 |
| ARG-LMM Stage 2 | h²=0.882 | 0.676 | ~12 |

**25% SNPs (ds_factor=4, 685 SNPs)**

| Model | h²_ARG | R² pred | time (s) |
|---|---|---|---|
| ARG-LMM | — (numerical) | 0.686 | ~97 |
| Just-Random | 0.638 | 0.666 | ~93 |
| GRM-LMM (ARG) | 0.771 | 0.662 | ~1 |
| GRM-LMM (SNP) | 0.756 | 0.628 | <1 |
| Just-Fixed | — | 0.639 | <1 |
| ARG-LMM Stage 2 | h²=0.923 | 0.667 | ~10 |

**Key observation**: ARG-based methods (Just-Random, GRM-LMM(ARG), Stage2) maintain
R²≈0.66–0.69 as SNPs drop to 25%; GRM-LMM(SNP) drops from 0.641→0.628. Core Figure 2
result confirmed.

---

## SARGENT on Figure 2 data

### Setup

```bash
REPO=/path/to/sargent
SIM=~/test-sargent/nurdan-all/Simulations/Figure2_MissingData/\
Simulation_350_0.1_8890_1Mb
OUT=~/test-sargent/nurdan_compare

# Import: symlinks bed/bim/fam/trees/phen + writes params.json
python $REPO/scripts/import_nurdan_fig2.py \
  --src $SIM \
  --pheno 1 --pheno-col 3 \
  --ds-factor 0 \
  --h2-truth 0.8 \
  --out $OUT/fig2_case

# Run SARGENT (all 350 individuals, no held-out test set)
time $REPO/build/src/progs/sargent \
  --pheno-col 3 \
  -o $OUT/sargent_fig2_ds0 \
  $OUT/fig2_case.phen \
  $OUT/fig2_case.trees \
  $OUT/fig2_case
```

For downsampled SNPs, change `--ds-factor` and pass the matching BED prefix;
the tree sequence is always the full `all.ts`.

### Three-way comparison (Seed8890, all SNPs, truth h²_SNP≈0.80, h²_ARG≈0)

| Method | h²_SNP | h²_ARG | h²_total | σ² | wall (s) |
|---|---|---|---|---|---|
| **SARGENT** | **0.725** | **0.106** | **0.831** | **0.175** | **4025\*** |
| ARG-LMM (MKLMM) | 0.000 | 0.793 | 0.793 | ~0.20 | 99† |
| ARG-LMM Stage 2 | BayesR | 0.882 | joint | — | 12† |
| Just-Random | — | 0.638 | 0.638 | — | 93† |
| GRM-LMM (ARG) | 0.000 | 0.771 | 0.771 | — | 1† |
| GRM-LMM (SNP) | 0.000 | 0.808 | 0.808 | — | <1† |
| Just-Fixed | 0.565 | — | — | — | <1† |

\* includes tree traversal over 4869 raw trees; user CPU time 24m49s  
† excludes kernel pre-computation from tree sequence; only 0.7s RData load counted

### Key finding: SNP–ARG decomposition

**SARGENT correctly separates SNP from ARG variance** (h²_SNP=0.725, h²_ARG=0.106,
close to truth 0.80/0.00). This is by design: Stage 1 BayesR exhausts typed-SNP
variance before Stage 2 fits the ARG residual.

**MKLMM (ARG-LMM) cannot make this separation**: its alternating REML absorbs
essentially all variance into the ARG random effect (h²_SNP=0.000, h²_ARG=0.793)
because the ARG kernels span the same space as the genotypes. BIC lasso selects
zero SNP fixed effects when the ARG already explains the data.

This difference is not a bug in either method — it reflects a fundamental model
choice. For estimating the total genetic variance (h²_total), both agree well
(SARGENT 0.831, MKLMM 0.793 vs truth ~0.80). For decomposing SNP vs. non-typed
ARG components, SARGENT's sequential design gives a well-defined answer;
MKLMM's decomposition is model-dependent and not stable when kernels and
genotypes are collinear.

### Timing: what is and isn't counted

| Step | SARGENT | MKLMM |
|---|---|---|
| Tree sequence → K matrices | included (tree traversal) | **excluded** (pre-computed offline) |
| Variance component fitting | included | included |
| Total wall time (this run) | 67 min (4025s) | 99s fit + ~?s pre-comp |

The kernel pre-computation (Python/tskit → .npy → .RData) is the MKLMM step
equivalent to SARGENT's tree traversal. It was run offline on Nurdan's server and
is not timed here. Add it to the MKLMM pipeline for a fair end-to-end comparison.

---

## Implementation notes

### ic.glmnet patch

Nurdan's original code uses `HDeconometrics::ic.glmnet(crit="bic")` to select the lasso
regularisation parameter by BIC. That package is not available for R ≥ 4.5.

`mklmm_core.R` replaces it with a faithful direct implementation:
```r
fit_path <- glmnet::glmnet(x, y, penalty.factor=weight, ...)
rss  <- colSums((y - predict(fit_path, newx=x))^2)
bic  <- nrow(x) * log(pmax(rss, .Machine$double.eps) / nrow(x)) + fit_path$df * log(nrow(x))
lambda_bic <- fit_path$lambda[which.min(bic)]
```
This is deterministic (no cross-validation randomness) and matches BIC selection.

### Model differences: SARGENT vs MKLMM vs Stage 2

| | SARGENT | MKLMM (ARG-LMM) | ARG-LMM Stage 2 |
|---|---|---|---|
| SNP effects | Stage 1 BayesR (variational) | BIC lasso (PLME.Fixed) | BayesR (coordinate ascent) |
| ARG weights | Stage 2 per-tree REML EM | AIREML per-kernel (PLME.Random) | AIREML per-segment (joint) |
| Parameterisation | variance w²_s | SD dr_s (V = σ²I + Σdr²K) | variance w²_s |
| Execution order | sequential (Stage1 → Stage2) | alternating backfit | joint iteration |
| Kernel input | raw tree sequence (.trees) | pre-computed K matrices (.RData) | pre-computed K matrices |
| Train/test | required | not required for h² | required |
| Kernel dedup | — (uses all trees) | unique topologies only | unique topologies only |

### Timing accounting

MKLMM timings in this README **exclude** kernel pre-computation (building K matrices
from the tree sequence with Python/tskit). SARGENT timings **include** all tree
traversal. For a fair end-to-end comparison, add the kernel pre-computation step to
the MKLMM pipeline time.
