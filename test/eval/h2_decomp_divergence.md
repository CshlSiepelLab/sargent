# Stage-1 `h²_snp` divergence on extreme cross-ancestry cells (diagnosis)

Date: 2026-05-31. Status: **diagnosed, not fixed** (production C — audit only,
per the guardrails). Surfaced by Block 1.1 v1 (PR #21). Tracked as a GitHub issue.

## Symptom

On 4 extreme cross-ancestry cells — `mix0.20_typed0.70` s1/s4/s5 and
`mix0.05_typed0.10` s4 — sargent fits a **nonsensical fixed-vs-random `h²`
decomposition**:

| cell | sim realised h²_snp / h²_arg | sargent FIT h²_snp / h²_arg |
|---|---|---|
| mix0.20_typed0.70_s1 | 0.297 / 0.297 | **−3.13 / +3.86** |
| mix0.20_typed0.70_s4 | 0.316 / 0.316 | **−3.03 / +3.03** |
| mix0.20_typed0.70_s5 | 0.289 / 0.289 | **−3.39 / +4.12** |
| mix0.05_typed0.10_s4 | 0.327 / 0.327 | **−1.16 / +1.80** |

The **simulated data is correct** (realised h² ≈ targets 0.30/0.30); the fit is
degenerate. Consequence: the BayesR fixed-effect betas are garbage, so the
BayesR and SARGENT point predictions collapse to ~0 `r2_calib` on these cells,
while ARG-BLUP-uniform (no fixed effect) is unaffected.

## Mechanism (Stage-1-only trajectory)

`sargent --snp-only` (Stage 1 = BayesR fixed-effect EM only) on
`mix0.20_typed0.70_s1`:

```
var_fixed iter  1: sigma2=0.2226  h2_snp=+0.7412  pi=[0.299 0.294 0.259 0.148]
var_fixed iter 50: sigma2=3.5560  h2_snp=-3.1355  pi=[0.527 0.299 0.009 0.165]
var_fixed iter 52: sigma2=3.5544  h2_snp=-3.1336  pi=[0.529 0.298 0.009 0.164]  (converged, rel_dsig 7e-5)
```

`sigma2` (residual variance) **inflates from 0.22 → 3.55**, i.e. to ~3.5× `var_y`
(≈1 for the standardised phenotype). Since `h2_snp = 1 − sigma2/var_y`
(`variational.c:361`), `sigma2 > var_y` ⇒ `h2_snp < 0`. It is a **stable
degenerate fixed point**, not a transient blow-up — the EM converges to it.
The mixture simultaneously loads ~16% of SNPs into the largest-effect class
(`pi[3]=0.164`) while collapsing the small-slab class (`pi[2]: 0.26→0.009`).

The relevant invariant is **violated**: `variational.c:347` asserts the
`σ² = ‖r‖²/n` update is *"bounded by var(y)"* (and only half-steps upward moves
to damp the known "mixture-snap" cascade — the comment notes `n=1000, p=3340`
diverged without it). At `n=500, p≈3100` in the mixed-ancestry regime the
half-step damping is **insufficient**: `‖r‖²/n` exceeds `var_y` (the fitted
betas overshoot — plausibly because the 20% pop-B training individuals, under
the mixed-population standardisation of allele-frequency-divergent SNPs, inflate
the residuals and feed the cascade), and nothing clamps `σ² ≤ var_y`.

## Why this is independent of the CG-cap fix (PR #21)

Stage 1 (`var_fixed`) is untouched by PR #21 (the cap is in `var_random`,
Stage 2). The uncapped run reproduces `h2_snp = −3.13` to the digit. So this is
pre-existing; the cap merely let the 4 cells finish (instead of hanging in
Stage 2) and thereby exposed it. It is the fixed-vs-random non-identifiability
(sargent's own `--help`: "the classical fixed-vs-random decomposition is not
identifiable") destabilising under cross-ancestry overfitting.

## Tested candidate fixes — all three obvious ones FAIL

Stage-1-only A/B tests on `mix0.20_typed0.70_s1` (truth h²_snp ≈ 0.30, baseline
fit −3.13):

| fix | resulting h²_snp | verdict |
|---|---|---|
| clamp σ² ≤ var_y | +0.000 | **counterproductive** — pegs reported h² at 0 but makes betas **14× larger** (‖β‖ 13 → 186). Holding σ² low keeps `u_bar = n/σ²` high ⇒ SNPs look *more* significant ⇒ *more* overfit. |
| upward-damping factor (0.5 → 1.0) | −3.13 | no effect |
| joint-EM (`--em`) | −3.13 | first `var_fixed` pass diverges identically; random effect inherits it |

**Conclusion: `h²_snp = −3` is a symptom; the disease is cross-ancestry
overfitting.** The fitted betas escape the prior entirely — BayesR largest
mixture class is σ² = 1e-2·var_y ⇒ prior std ≈ 0.093, yet the degenerate cell
fits **max|β| = 11.9** (~128 prior-stds; **40×** the well-behaved cells'
max|β| = 0.31). On mixed-ancestry training (80% A / 20% B) the coordinate
descent fits allele-frequency-divergent SNPs to separate ancestry
(confounding-by-ancestry); the EM then raises σ² to ~4× var_y as a partial
shrinkage response. The high baseline σ² is the EM's own (partial) correction —
which is why clamping it *down* makes things worse.

## Recommended direction (methodological — for Adam, none applied)

A guard/clamp cannot fix this; it needs a real intervention, all of which
change behaviour on **all** cells and so are deferred for decision:

1. **Ancestry-PC covariates** in Stage 1 — standard PRS confounding control.
2. **Cap total fixed-effect variance / sparser-or-capped large-class prior** so
   betas cannot escape the prior.
3. **Constrain the coordinate descent** (bound per-SNP β, or proper p ≫ n
   shrinkage).
4. **Filter / down-weight ancestry-divergent SNPs.**

For Block 1.1 the 4 affected cells are excluded from the verdict and the
conclusion is unchanged. Tracked as GitHub issue #22.

## Reproduce

```bash
D=$HOME/test-sargent/block1_1_v1/mix0.20_typed0.70_s1
build/src/progs/sargent --snp-only --nprobe 10 --outroot /tmp/d/est \
  "$D/case.phen" "$D/case.trees" "$D/case"   # watch var_fixed sigma2 climb past var_y
```
