# Native out-of-sample prediction in the `sargent` binary (scope 3b)

**Status:** in progress on branch `adam/native-predict`. Implements the "full"
scope (3b) of the proposal `native_prediction_proposal.md` (in the paper repo):
`sargent` itself computes the out-of-sample predictions and writes
`<outroot>.predict.json`, so the paper's numbers come from `sargent` rather than
from a Python reimplementation.

**Decision (Adam, 2026-07-30):** full ladder parity — reproduce all five arms of
`analysis/predict.py` (`bayesR`, `snp_lmm`, `arg_lmm`, `sargent_u`, `sargent`)
so `est.predict.json` is reproducible and `predict.py` can eventually be retired.

**Merge gate (Adam):** must match the current `predict.py` path to < 1e-6 on
every metric before merge. See `test/eval/native_predict_test.sh`.

---

## 1. What `sargent --predict` does

After the normal fit on the **train** ARG (`case.trees`, `case.bed`,
`case.phen`) — Stage 1 gives β̂, Stage 2 gives ŵ²[s] — the prediction stage:

1. Loads the **full** ARG `case.full.trees` (train+test, individuals
   `[0,n_train)` = train, `[n_train,n_total)` = test).
2. Maps ŵ²[s] from train-tree index → full-tree index by midpoint
   (`map_train_to_full`; no pruning by default, so `sg->ws2[s]` = train tree s).
3. Builds four ARG kinship blocks over the full ARG in **one** set of `n_train`
   tree passes (unit-vector columns over train individuals; each pass yields a
   full-length `G_s·e_j`, whose train rows feed the train block and test rows the
   cross block), accumulating **two** weightings simultaneously:
   - uniform (span-weighted eGRM): `K_uni = Σ_s span_w[s]·G_s`
   - learned:                      `K_sg  = Σ_s ŵ²[s]·span_w[s]·G_s`
     where `span_w[s] = span_s / mean_span` (the walker's `st.span`,
     `use_span=1`), identical to `predict.py`'s `span_w`.
   Blocks: `K_*_train` (n_train×n_train), `K_*_cross` (n_test×n_train).
4. Loads genotypes and standardizes **exactly as predict.py** (empirical TRAIN
   mean and population std, ddof=0), aligned to sargent's train-polymorphic
   variant set (β-indexed): `X_train`, `X_test` (shape p×n).
5. Typed-SNP GRM: `K_g_train = Xtrainᵀ Xtrain / p`, `K_g_cross = Xtestᵀ Xtrain / p`.
6. Applies the five predictors (below), computes metrics, writes JSON.

Genotypes come from BED (as in predict.py). `--ts-geno` + `--predict` is rejected.

## 2. Exact math to reproduce (from `analysis/predict.py`)

**center_kinship(K_train, K_cross):**
```
r = rowmean(K_train)        # (n_tr,1)
c = colmean(K_train)        # (1,n_tr)
g = grandmean(K_train)      # scalar
K_train_c = K_train - r - c + g
d = mean(diag(K_train_c));  if d>0: K_train_c /= d
# cross uses TRAIN col/grand means, its OWN row means:
K_cross_c = K_cross - rowmean(K_cross) - c + g;  if d>0: K_cross_c /= d
```
Double-centering annihilates the rank-1 constant offset in the ARG kinship;
mean-diag normalization fixes global scale. Both blocks share `c`, `g`, `d` so
they stay on the same scale (required for BLUP). Global scale/additive offsets
in our kernel vs predict.py's cancel here → C blocks need only match up to those.

**fit_blup(y, K):** single-random-effect REML EM, `h2_init=0.5`, `max_iter=50`,
`tol=1e-5`.
```
y_c = y - mean(y);  var_y = pop_var(y_c)
su2 = 0.5*var_y;  s2 = 0.5*var_y
repeat up to 50:
  V = su2*K + s2*I
  L = chol(V)            # on failure: V += 1e-6*I; chol again
  Vinv_y = V^{-1} y_c    # via L
  Vinv   = V^{-1}        # full inverse via L
  tr_VinvK = Σ_ij Vinv⊙K ; tr_Vinv = trace(Vinv)
  a_u = (K Vinv_y)·Vinv_y ; a_e = Vinv_y·Vinv_y
  new_su2 = su2 * a_u / max(tr_VinvK,1e-30)   # EM multiplicative
  new_s2  = s2  * a_e / max(tr_Vinv, 1e-30)
  clamp each to >= 1e-8*var_y
  if both rel-changes < tol: break
lam = s2 / max(su2,1e-30)
alpha = (K + lam·I)^{-1} y_c
return s2, su2, alpha, mean(y)
```

**Predictors** (mean_y = train mean; r = y_train − Xtrainᵀβ):
- `bayesR`   : ŷ = Xtestᵀβ + mean_y
- `snp_lmm`  : center K_g; fit_blup(y_train,K_g_c); ŷ = K_g_cross_c·α + mean_y;
              `snp_lmm_h2 = su2/(su2+s2)`
- `arg_lmm`  : center K_uni; fit_blup(y_train,K_uni_c); ŷ = K_uni_cross_c·α + mean_y;
              `arg_lmm_h2`
- `sargent_u`: r; center K_uni; fit_blup(r,K_uni_c)→α,mean_r; ŷ = Xtestᵀβ + mean_r + K_uni_cross_c·α
- `sargent`  : same with K_sg; also `sargent_sigma2`; if `case.test.true_components.tsv`
              present, `sargent_rand_r = corr(ŷ − Xtestᵀβ, random_true)`

**metrics(ŷ,y):** `r=corr`, `r2=1−SSres/SStot`, `r2_calib=r²`, `mse=SSres/n`.

**JSON keys** (match predict.py exactly): `case, n_train, n_test,
n_trees_full, n_trees_train, bayesR, snp_lmm, snp_lmm_h2, arg_lmm, arg_lmm_h2,
sargent_u, sargent, sargent_sigma2, sargent_rand_r`. Metric objects:
`{r, r2, r2_calib, mse}`.

## 3. CLI

```
sargent <train.phen> <train.trees> <geno_prefix> \
        --predict \
        [--full-trees PATH]   # default <geno_prefix>.full.trees
        [--test-geno PREFIX]  # default <geno_prefix>.test
        [--test-pheno PATH]   # default <geno_prefix>.test.phen
        [--predict-out PATH]  # default <outroot>.predict.json
```
`--predict` runs Stage 1 + Stage 2 as usual (needs β̂ and ŵ²), then the
prediction stage. `n_test` is inferred as `num_individuals(full.trees) − n_train`.
`n_test == 0` (no test set) is a hard error under `--predict`.

## 4. Files touched

- `include/geno.h`, `src/lib/bed.c` — expose per-poly-variant file index and a
  raw-read-by-file-index helper for aligned test-BED reads.
- `include/predict.h`, `src/lib/predict.c` — the whole prediction pipeline.
- `src/progs/sargent.c`, `src/progs/sargent.help.h`, `src/lib/CMakeLists.txt`.

## Result (merge gate — PASS)

`native_predict_test.sh` on a small panmictic case (n_train=150, n_test=60,
L=300 kb, 679/717 train/full trees): **all 29 shared metrics agree with
`predict.py` within 1e-6** (worst relative diff 4.8e-8 at `sargent.r2_calib`).

- `arg_lmm`, `snp_lmm` agree to ~1e-13 (C branch kernel == tskit/GRM path after
  `center_kinship`; deterministic linear algebra).
- `sargent`, `sargent_u`, `bayesR` agree to ~1e-8 — floating-point summation
  order (Python streaming loop vs C incremental walker) and REML-EM iteration
  noise; well under the gate.

T5 (no-op without `--predict`): `.beta.tsv`/`.w2.tsv` byte-identical, no
`.predict.json` written.

## Reproduce (merge gate)

```bash
# from ~/git-src/sargent, sargent env active, binary built
bash test/eval/native_predict_test.sh
# generates a small case, runs sargent --predict and predict.py,
# asserts every est.predict.json metric agrees < 1e-6, and that n_test=0 is a no-op.
```
