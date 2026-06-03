# SARGENT — overview for Nurdan

This document explains how to run SARGENT on your simulated data, the
major operating modes, and how it relates to MKLMM / ARG-LMM.  For
full CLI reference run `sargent --help`.  For side-by-side comparison
results see `README.md` in this directory.

---

## 1. Build

```bash
git clone --recurse-submodules https://github.com/asiepel/sargent.git
cd sargent
mamba activate sargent     # provides cmake, meson, ninja, tskit, msprime
cmake -B build -S .
cmake --build build -j
./build/src/progs/sargent --help
```

---

## 2. Converting your simulation directories to SARGENT format

Two thin adapter scripts create a `case.*` layout from your directory
structure — no data is copied, just symlinks + a translated
`params.json`.

### Figure 1 simulations (all N individuals, no train/test)

```bash
python scripts/import_nurdan.py \
    --src ~/test-sargent/nurdan-all/Simulations/Figure1_Simulation_Xb_U/\
Simulation_Fig1_200ind_CausalRatio_0.1_Seed6903_1Mb \
    --trait vc_balanced \
    --pheno-col 1 \
    --out /tmp/case
```

This writes `case.{phen,trees,bed,bim,fam,params.json}`.  The adapter
prints the exact `sargent` invocation to use, including the right
`--pheno-col` offset (your `.pheno` files have two ID columns before the
replicate columns, so phenotype replicate 1 is at column 3).

### Figure 2 simulations (SNP-downsampled, 350 individuals)

```bash
python scripts/import_nurdan_fig2.py \
    --src ~/test-sargent/nurdan-all/Simulations/Figure2_MissingData/\
Simulation_350_0.1_8890_1Mb \
    --pheno 1 --pheno-col 3 \
    --ds-factor 0 \          # 0=100% SNPs; 4≈25%; 10≈10%
    --h2-truth 0.8 \
    --out /tmp/case_fig2
```

---

## 3. Basic run

```bash
sargent \
    --pheno-col 3 \
    --outroot /tmp/myrun \
    /tmp/case.phen \
    /tmp/case.trees \
    /tmp/case          # BED prefix: case.bed / case.bim / case.fam
```

By default SARGENT runs the **sequential** pipeline (Stage 1 → Stage 2;
see §4) and writes three output files (see §5).  Progress is printed to
stderr; redirect with `2>run.log` or suppress with `--silent`.

A `.gall` cache file is written next to the BED prefix on the first run
(eigendecomposition of the aggregate ARG kinship); it is reused
automatically on subsequent runs of the same dataset.

---

## 4. Operating modes

| flag | what runs | use when |
|---|---|---|
| *(none — default)* | Stage 1 (BayesR) → Stage 2 (per-tree REML EM) | standard fit |
| `--snp-only` / `-S` | Stage 1 only | check typed-SNP baseline; fast |
| `--arg-only` / `-A` | Stage 2 only (β = 0) | ARG-only fit, no fixed effects |
| `--uniform-arg` / `-U` | Stage 2 with uniform weights (single τ²) | ARG-BLUP-uniform; fastest Stage 2 |
| `--em N` | sequential then N joint-EM polish iterations | squeeze a few % more h²_arg on small-n cases |

### Combinations that make sense for your comparisons

```bash
# Equivalent to MKLMM "Just-Fixed" (lasso on SNPs, no ARG):
sargent --snp-only ...

# Equivalent to MKLMM "Just-Random" / "GRM-LMM(ARG)" (ARG only, uniform):
sargent --arg-only --uniform-arg ...

# Closest to MKLMM "ARG-LMM" (ARG with per-kernel weights, no typed SNPs):
sargent --arg-only ...           # uses SARGENT's per-tree REML EM

# Full SARGENT hybrid (BayesR SNPs + per-tree-weighted ARG):
sargent ...                      # default
```

### Notes on `--em N`

Joint EM warm-starts from the sequential posterior and runs N additional
iterations of interleaved (β, w²_s, σ²) updates.  It is useful on
small-n cases (n ~ 200) where it modestly improves the ARG–SNP split.
At larger n or with many trees the improvement is minimal and the
divergence-revert safety net usually fires.  **Use `--em 10` or `--em 20`
on your data; do not use `--em 50`** — the default 50-outer-iteration
cap is too aggressive for cases with thousands of trees.

---

## 5. Output files

| file | contents |
|---|---|
| `<outroot>.summary.json` | h²_snp, h²_arg, h²_total, σ², n_indiv, n_variants, n_trees |
| `<outroot>.beta.tsv` | per-SNP posteriors: mixture responsibilities γ_jk and E[β_j] |
| `<outroot>.w2.tsv` | per-tree posteriors: responsibilities ρ_sm and E[w²_s] |

**`summary.json` is the authoritative variance-component output.**
Parse it with `json.load(open(...))` in Python or `jsonlite::read_json()`
in R.  The `.beta.tsv` and `.w2.tsv` files have one row per variant /
tree respectively.

The decomposition is defined sequentially:
- `h2_snp` = fraction of `var(y)` explained by the Stage-1 BayesR fit
  on typed SNPs (= 1 − var(residual)/var(y) after Stage 1).
- `h2_arg` = fraction of `var(y)` additionally explained by the
  Stage-2 ARG random effect on the residual.
- `h2` = `h2_snp + h2_arg`.

---

## 6. How SARGENT compares to MKLMM / ARG-LMM

### What is the same

- **Model**: y = Xβ + u + ε, u ~ N(0, Σ_s w²_s G_s), ε ~ N(0, σ²I).
  SARGENT and MKLMM both estimate per-tree (per-kernel) weights w²_s
  and variance components.
- **Inputs**: both read a tskit `.trees` file and PLINK BED genotypes.
- **Uniform-weight special case**: SARGENT's `--uniform-arg` mode is
  numerically equivalent to MKLMM's "GRM-LMM(ARG)" (single averaged
  ARG kernel, EMMA/AIREML solve).
- **Heritability ballpark**: on matched inputs total h² agrees to
  within ~5% (see comparison table in `README.md`).

### What is different

| aspect | SARGENT | MKLMM (ARG-LMM) |
|---|---|---|
| **K matrix handling** | traverses the tree sequence at runtime per-iteration; no pre-computation | K matrices pre-computed offline once, reused cheaply |
| **SNP fixed effects** | BayesR (variational spike-and-slab; 4 mixture classes) | BIC lasso (glmnet path + BIC selection) |
| **Weight inference** | multiplicative REML EM per tree (all S trees retained) | AIREML per unique-topology kernel (deduplicates to ≪ S kernels) |
| **Execution order** | sequential: Stage 1 exhausts SNP variance, then Stage 2 fits the residual | alternating backfit on (β, w²_s) jointly |
| **SNP/ARG split** | well-defined by construction (sequential priority) | not separable when SNP and ARG kernels are collinear (model-choice dependent) |
| **Train/test split** | required (SARGENT computes OOS predictions) | not required for h² estimation |
| **Timing** | includes all tree traversal | usually quoted excluding K pre-computation — add that for a fair comparison |
| **Language** | C binary | R scripts |

### The decomposition difference in practice

On a pure-SNP trait (truth h²_SNP = 0.8, h²_ARG = 0):

| | h²_SNP | h²_ARG | h²_total |
|---|---|---|---|
| SARGENT (sequential) | 0.73 | 0.11 | 0.83 |
| MKLMM (ARG-LMM) | 0.00 | 0.79 | 0.79 |
| truth | 0.80 | 0.00 | 0.80 |

Neither number is "wrong" — they reflect different model conventions.
SARGENT's sequential design gives a stable SNP-prioritised answer;
MKLMM's alternating fit absorbs variance into whichever component fits
first.  For total h² and prediction R² both methods agree well.  See
`README.md` for the full three-way (Seed8890) comparison.

---

## 7. Tips and gotchas

**Phenotype column.** Your `.pheno` files have two ID columns then K
replicate phenotypes.  Replicate j is at `--pheno-col (j + 2)`.
Example: replicate 1 → `--pheno-col 3`.

**Large runs / server.** Run in a tmux session so it survives SSH drops.
SARGENT is resumable if you keep the `.gall` cache in place (no need to
re-compute the eigendecomposition).

**Cross-ancestry cells.** If you run on mixed-ancestry training cohorts,
the Stage-1 BayesR can overfit allele-frequency-divergent SNPs and
produce `h2_snp < 0` on extreme cells.  The binary has a trust-region
cap (`SARGENT_S2_VFCAP=0.6`, enabled by default) that recovers the h²
decomposition; it does not recover prediction on those cells.

**Unique-topology deduplication.** MKLMM deduplicates trees to unique
topologies before building K matrices, which greatly reduces the
effective number of kernels.  SARGENT operates on all S trees in the
sequence.  The practical effect is that SARGENT's per-tree weights are
noisier per-segment but cover the full sequence without any topology
binning.

**Fused-lasso flags (`--ws2-l1`, `--ws2-l2`).** Both default to 0
(legacy behaviour, byte-identical to the unregularised per-tree EM).
Setting `--ws2-l1 0.5` modestly improves **localization AUC** in
oligogenic regimes (active trees form a few concentrated blocks) but
does not improve heritability recovery or prediction.  Leave at 0 for
standard variance-component analyses.
