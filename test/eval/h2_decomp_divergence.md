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

## Candidate fixes (for Adam — do not modify production C without sign-off)

1. **Hard floor `h2_snp ≥ 0` / clamp `σ² ≤ var_y`** — enforces the invariant the
   comment already assumes. Cheap and makes the *reported* h² sane, but does
   NOT fix the underlying overfit betas, so cross-ancestry prediction stays bad.
2. **Stronger σ² cascade control** in `var_fixed` (symmetric damping / clip, as
   `var_joint` uses), or a divergence guard that detects `σ² → >var_y` and
   backtracks. Targets the root EM instability.
3. **Joint-EM mode (`--em N`)** for these cells — let the ARG random effect
   absorb the h²_arg variance so Stage 1 stops over-attributing it to (and
   overfitting) the typed SNPs. Worth checking whether `--em` already avoids the
   divergence on these cells (quick test, no code change).
4. **Investigate the mixed-ancestry standardisation / mixture-assignment** — why
   16% of SNPs land in the top effect class on mixed training. This is the
   cross-ancestry analogue of the BayesR fixed-effect blow-up that also drives
   the negative raw R² (see the `r2_calib` patch rationale in PR #21).

## Reproduce

```bash
D=$HOME/test-sargent/block1_1_v1/mix0.20_typed0.70_s1
build/src/progs/sargent --snp-only --nprobe 10 --outroot /tmp/d/est \
  "$D/case.phen" "$D/case.trees" "$D/case"   # watch var_fixed sigma2 climb past var_y
```
