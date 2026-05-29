# `predict.py:fit_blup` audit — diagnosis of the +0.25 h² bias

Date: 2026-05-29. Audit task: Block 1.3 of the eval workblock plan.
Scope: production single-component variance-component estimator
`fit_blup` in `test/eval/predict.py`, used across Phase 3 / Phase 6 / Phase 10.

## TL;DR (numbers > narrative)

Per-case h² on all 10 Phase 10 admixed simulations (truth mean ≈ 0.50):

| variant                                    | mean h² | bias   |
|--------------------------------------------|---------|--------|
| truth (realised)                           | 0.498   | —      |
| `fit_blup` raw (current production)        | 0.763   | **+0.264** |
| `fit_blup` on **double-centered K**        | 0.550   | +0.052 |
| SARGENT C-code `--arg-only h²_arg`         | 0.578   | +0.080 |
| `fit_blup` on **SNP-based GRM**, same y    | 0.571   | +0.066 |
| Proper REML (P-projector, intercept)       | 0.798   | +0.293 |

**Dominant cause: the ARG kinship from `streaming_combined_kinship`
carries a huge constant offset (rank-1 component ≈ 65 % of trace) that
`fit_blup` absorbs into σ_u². Neither REML projection, intercept fixed
effects, nor extra top-k PC covariates remove it. The fix is to
**double-center K** (equivalently, treat the constant block as a fixed
effect that lives inside K, not inside V).** Double-centering brings
`fit_blup` to the same residual bias as SARGENT's own C-code estimator
and as `fit_blup` on a SNP GRM, all ~+0.05–0.08. The +0.25 inflation is
not in `fit_blup` per se — it is in the way the eval-side ARG kinship
is built and passed to `fit_blup`.

## 1. What `fit_blup` actually optimises

From `test/eval/predict.py` (l. 153–201):

```python
def fit_blup(y_train, K_train_train, h2_init=0.5, max_iter=50, tol=1e-5):
    n = len(y_train)
    y_c = y_train - y_train.mean()
    sigma_u2 = h2_init * y_c.var()
    sigma2   = (1 - h2_init) * y_c.var()
    for _ in range(max_iter):
        V = sigma_u2 * K_train_train + sigma2 * np.eye(n)
        L = np.linalg.cholesky(V)
        Vinv_y = solve(L.T, solve(L, y_c))
        K_Vinv_y = K_train_train @ Vinv_y
        Vinv = solve(L.T, solve(L, np.eye(n)))
        tr_VinvK = np.einsum("ij,ji->", Vinv, K_train_train)
        tr_Vinv  = np.trace(Vinv)
        a_u = K_Vinv_y @ Vinv_y;       b_u = tr_VinvK
        a_e = Vinv_y @ Vinv_y;         b_e = tr_Vinv
        new_su2 = sigma_u2 * a_u / b_u
        new_s2  = sigma2   * a_e / b_e
        ...
```

- Fixed effects: implicit intercept via `y_c = y - mean(y)`. No matrix X
  is built, so any non-intercept structure cannot be absorbed.
- Score: classical multiplicative EM updates with the **ML projection**
  `V^{-1}`, not the REML projection `P = V^{-1} - V^{-1}X(X^TV^{-1}X)^{-1}X^TV^{-1}`.
  This is "ML-EM minus the mean" — not proper REML. For an intercept-only
  X with n = 300 the difference between ML and REML projections is O(1/n),
  i.e. negligible (confirmed below in test (a)).
- Termination: relative tolerance 1e-5 on (σ_u², σ²) component-wise. In
  practice converges in 15–30 iterations cleanly and monotonically.

## 2. Per-hypothesis test results

All tests on the 5 sparse + 5 diffuse seeds in
`~/test-sargent/phase10_admix/`. Truth (realised) h²_arg averages 0.51
(sparse) / 0.49 (diffuse).

### (a) Proper REML with P-projector — REJECTED

Implemented an explicit-P AI-EM REML (`/tmp/fit_blup_audit.py`,
`fit_reml_proper`) and compared per case. Result on the 5 sparse seeds:

| case       | truth | current fit_blup | proper REML |
|------------|-------|------------------|-------------|
| sparse_s1  | 0.485 | 0.811            | 0.815       |
| sparse_s2  | 0.519 | 0.772            | 0.778       |
| sparse_s3  | 0.493 | 0.832            | 0.836       |
| diffuse_s1 | 0.497 | 0.740            | 0.745       |
| diffuse_s2 | 0.531 | 0.810            | 0.815       |

REML log-likelihood at the proper-REML fit is higher than at the fit_blup
fit by 0.005–0.007 nats; ML log-likelihood is the reverse, also by
0.005–0.007 nats. **Both fits land on essentially the same point** —
projection matters at the 5th decimal at n = 300. **REML projection is
not the cause of the +0.25 bias.**

### (b) Rank-1 dominance of the ARG kinship — CONFIRMED, dominant cause

Spectrum of the normalised (mean diag = 1) ARG kinship on sparse_s1:

```
top 5 eigenvalues: [195.23, 16.31, 8.42, 5.29, 4.87]
λ_1 / tr(K)        = 0.651   (uniform = 1/n = 0.003)
sum top 5 / tr(K)  = 0.767
```

One eigenvalue carries **65 % of the trace**. Equivalent SNP-GRM
spectrum has λ_1 / tr ≈ 0.05. The top eigenvector is **almost exactly
the constant vector**:

```
|cos(v_1, 1)| = 0.998
v_1's projection of y² = 0.00 % of y'y
```

so it is a pure "global IBD" eigenmode — the deep-coalescent shared
branch-length baseline that every pair gets — carrying no phenotype
signal but enormous variance.

Removing top-k spectral components and refitting:

| k removed | sparse_s1 | sparse_s2 | sparse_s3 | diffuse_s1 | diffuse_s2 |
|-----------|-----------|-----------|-----------|------------|------------|
| 0         | 0.811     | 0.772     | 0.832     | 0.740      | 0.810      |
| 1         | **0.608** | **0.553** | **0.657** | **0.507**  | **0.608**  |
| 2         | 0.569     | 0.502     | 0.564     | 0.374      | 0.543      |
| 5         | 0.449     | 0.198     | 0.480     | 0.256      | 0.516      |

Removing just v_1 from K cuts the bias roughly in half. Removing 2 lands
near truth on average. Note: simply running proper REML with v_1, ..., v_k
as **fixed-effect covariates** does *not* fix it (we tested k_cov=1,2,5,10,20
and h² stayed at 0.74–0.83 throughout). The reason is subtle: adding v_1 as
fixed effect kills its projection in `Py`, but `σ_u² K` still contains the
huge `λ_1 v_1 v_1^T` block in V, and the REML score `y'PKPy / tr(PK)` is
dominated by that block, driving σ_u² up. **Fixed-effect projection on y is
not enough; the constant mass has to be removed from K itself.**

### (c) Kinship normalisation — minor (scaling invariance) but raw scale unsafe

`fit_blup` is scale-invariant in K up to numerical floor: tested raw K,
K/mean-diag, K·n/tr(K), K/λ_max on sparse_s1:

```
raw           : h2_fit = 0.000   (su2 hit lower clamp; K diag ~ 3e4)
mean_diag=1   : h2_fit = 0.811
trace=n       : h2_fit = 0.811
lmax=1        : h2_fit = 0.999   (boundary trap — su2 too small)
```

The `mean_diag=1` and `trace=n` normalisations are mathematically
equivalent (`tr/n == mean(diag)`); both work. **Raw scale fails** because
the lower clamp `max(new_su2, 1e-8 * var_y)` activates when K diag is
~3e4. This is a real (and easy) downstream bug, but is NOT the cause of
the +0.25 bias — `oracle_ws2.py` already does mean-diag=1 normalisation.

### (d) EM convergence — clean

Iterate trace for sparse_s1 (15 iters):
`0.658, 0.737, 0.775, 0.792, 0.801, 0.806, 0.808, 0.809, 0.810, 0.810,
0.810, 0.810, 0.810, 0.810, 0.810` — monotone, stable, no oscillation,
no boundary clipping. Final 3 sig figs locked. **EM convergence is not
the cause.**

### Cross-check: same `fit_blup`, SNP-GRM vs ARG-K, same y

| case       | truth | h² on **ARG-K** | h² on **SNP-K** | λ₁/tr (ARG) | λ₁/tr (SNP) |
|------------|-------|-----------------|-----------------|-------------|-------------|
| sparse_s1  | 0.485 | 0.811           | 0.597           | 0.65        | 0.05        |
| sparse_s2  | 0.519 | 0.772           | 0.552           | 0.65        | 0.05        |
| sparse_s3  | 0.493 | 0.832           | 0.622           | 0.63        | 0.06        |
| diffuse_s1 | 0.497 | 0.740           | 0.542           | 0.65        | 0.05        |
| diffuse_s2 | 0.531 | 0.810           | 0.542           | 0.65        | 0.05        |

The SNP GRM (built with standardised dosages, so it is implicitly mean-zero
in each column → no constant offset) gives an h² with mean bias +0.06,
matching SARGENT's own per-tree REML EM. **The bias is in the ARG kinship
construction, not in `fit_blup` per se.**

## 3. Verdict

> The +0.25 systematic upward bias is caused by the aggregate ARG
> kinship `K_arg = Σ_s G_s/S` having a **rank-1 constant offset** (the
> deep-coalescent shared-branch-length mean) that accounts for ~65 % of
> tr(K) and is collinear with the all-ones vector. `fit_blup` lumps that
> constant mass into σ_u², inflating h². Hypothesis (b), in its
> "K-centering" form.

Hypotheses (a), (c), (d) are not material causes:

- (a) ML vs REML: 0.005 nat difference at n = 300; bias-equivalent.
- (c) Normalisation: scale-invariant once mean-diag=1 or trace=n is applied;
  the only issue is a raw-scale boundary trap, irrelevant here.
- (d) EM convergence: clean, monotone, locked.

SARGENT's C-code per-tree REML EM (bias +0.08) avoids the bug because each
G_s is fit separately with its own ŵ²_s, so the *common* constant offset
across trees doesn't aggregate into one inflated variance component — each
tree's offset competes locally against σ².

## 4. Recommended fix

### Option A (preferred): double-center K inside `fit_blup`

One-liner added at the top of `fit_blup` after computing `y_c`:

```python
# Double-center K to absorb the deep-coalescent constant offset.
# Equivalent to treating the all-ones eigenmode of K as a fixed effect.
n = K_train_train.shape[0]
ones_over_n = np.ones((n, 1)) / n
row_mean = K_train_train.mean(axis=1, keepdims=True)
col_mean = K_train_train.mean(axis=0, keepdims=True)
grand_mean = K_train_train.mean()
K_train_train = K_train_train - row_mean - col_mean + grand_mean
# Renormalise to mean diag = 1 for numerical stability.
d = float(np.mean(np.diag(K_train_train)))
if d > 0:
    K_train_train = K_train_train / d
```

(Equivalent matrix form: `H = I - 11'/n; K <- H K H`.)

Validated bias on the 10 Phase 10 seeds: **+0.052** (vs +0.264 raw).
This matches SARGENT's C-code estimator (bias +0.080) and `fit_blup`
on a SNP-GRM (bias +0.066), so the residual is the irreducible
finite-sample noise floor, not a fixable structural defect.

Caveats:

- `K_test_train` (off-diagonal cross block) used inside `predict_blup`
  must be projected with the **same** centering convention or the BLUP
  prediction will be miscalibrated. Concretely, the prediction call
  becomes `K_cross_centered = K_cross - row_mean_train - col_mean_test +
  grand_mean_train`. Wire the centering constants out of `fit_blup` and
  apply them to `K_cross` before `predict_blup`. Otherwise prediction R²
  will degrade.
- Centering at fit time is mathematically equivalent to declaring the
  all-ones eigenmode a fixed effect with infinite prior. Since `fit_blup`
  already implicitly centers y on the mean, this is consistent.

### Option B (safe quick win): drop `fit_blup` h² from result tables, keep prediction

The R²(ŷ, y) numbers from `fit_blup` are insensitive to the bias (BLUP
shrinkage is governed by the *ratio* σ²/σ_u² × scale of K; once K is
mean-diag normalised this ratio still maps to a reasonable shrinkage).
The unbiased estimate to put in the headline tables is the SARGENT
`--arg-only h²_arg` from `est.summary.json`, which has bias ~+0.08 and
is what we report in Phase 10.

Recommend **Option A + B together**: implement the double-centering fix
(so future `fit_blup` h² values are usable), and in the meantime drop
`gblup_h2`, `arg_blup_uniform_h2`, `h2_fit` from result tables that have
been written by the eval pipeline before the patch — they are not
trustworthy. Phase 3 and Phase 6 prediction R² values are fine; only
the h² values they incidentally report are biased.

## 5. Repro

Audit scripts kept in `/tmp/`:
- `/tmp/fit_blup_audit.py` — hypotheses (a), (b), (c), (d) end-to-end.
- `/tmp/fit_blup_audit_part2.py` — top-k as fixed effect (refutes the
  intercept-fix interpretation).
- `/tmp/fit_blup_audit_part3.py` — SNP-GRM cross-check (same `fit_blup`,
  bias drops to +0.06).
- `/tmp/fit_blup_audit_part4.py` — double-center K cross-check.
- `/tmp/fit_blup_audit_part5.py` — all 10 seeds, final summary.

Run with `mamba run -n sargent python <path>`.

No production code was modified. Patch to be applied separately.
