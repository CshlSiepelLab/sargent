# Block 2.6.5 — Spike-and-slab on `ws²` directly (design)

Date: 2026-06-01 (pm). Status: design; implementation pending.

## Goal

Implement the **mixture-prior posterior for `ws²_s`** that the SARGENT
`main.tex` Methods derives (eq. ~\ref{eq:rho-update}, around line 1241) —
the prior the code abandoned for point-estimate REML EM.

```
ws²_s | c_s=m ~ δ(τ_m²),    p(c_s = m) = φ_m,   Σ_m φ_m = 1
```

with `τ_1²=0` (spike) and `τ_2², ..., τ_M² > 0` (slab[s]). After fitting,
each tree gets a posterior responsibility `ρ_{s,m} ∈ [0, 1]` and
`E[ws²_s] = Σ_m ρ_{s,m} τ_m²`.

**Target outcome**: estimator M_e tracks truth M_e within 2× across the
Block 2.5.2 polygenicity sweep (currently flat at 12-20 regardless of
truth from 5 to 600).

**Critical caveat (Adam 2026-06-01)**: spike-and-slab might help
**localization** (concentrate ws² on truly active trees) but **hurt
prediction** (aggressive sparsification can throw away small-but-real
signal). The acceptance criteria below include a **no-regression-on-
prediction** test alongside the localization improvements.

## The challenge — already documented failure of the naïve approach

The original variational E-step (main.tex line 1241):

```
ρ_{s,m} ∝ φ_m · exp[(τ_m²/2)(A_s − B_s)]
```

where `A_s = r̃ᵀ G_s r̃`, `B_s = tr(V⁻¹ G_s)`, `r̃ = V⁻¹(y − Xβ̂)`.

**Documented failure (main.tex line 1253):** in the `S ≫ n` regime,
`(τ_m²/2)(A_s − B_s)` is of order `n²/S²` (≈ 10⁻³ at n=100, S≈10³),
which is too small to push the softmax responsibilities away from the
prior mean. `ρ` remains near-uniform across mixture components and σ²
fails to converge correctly. This is the failure mode that motivated the
point-estimate REML EM shortcut in the current code.

**For Block 2.6.5 to succeed where the original failed**, we need a
posterior approximation that gives a *stronger* per-segment score than
`τ_m²·(A_s − B_s)`. Three candidate paths:

### Option A — Wakefield-style Bayes factor for segments

The SNP-side BayesR update (main.tex line 1185) uses the Wakefield
approximation:

```
γ_{j,k} ∝ π_k · 𝒩(z_j; 0, u_j σ_k² + 1)
```

where `u_j = xⱼᵀ V⁻¹ xⱼ` is the whitened squared norm of SNP j.
This is *strongly discriminating* because `u_j` is O(n), so the
Bayes factor's variance ratio (`u_j σ_k² + 1`) is large enough to
separate active from null even at modest effect sizes.

The segment-side analogue would be a Bayes factor for segment s's
variance-component contribution:

```
BF_m(s) ≈ density ratio under N(0, τ_m² G_s) vs N(0, ε G_s)
ρ_{s,m} ∝ φ_m · BF_m(s)
```

Concretely, the segment-level score under a mixture component τ_m² is
approximately a noncentral χ² with d.o.f. ≈ rank(G_s) and non-centrality
proportional to A_s / B_s. The exact Bayes factor over the prior δ(τ_m²)
is computable in closed form via the variance-component LRT.

**Expected advantage**: per-segment score has scale O(rank(G_s)) ≈ O(n),
not O(n²/S²). Strong discrimination.

### Option B — Empirical-Bayes mixture-weight M-step

Update `φ_m` from data (full EM on the mixture proportions) instead of
holding it at the prior. As `φ_m` shifts toward the empirical data
distribution, the softmax pull becomes stronger.

Cheaper than (A); could be implemented as a 5-line addition. May not be
sufficient on its own — the underlying scoring issue remains — but
could combine well with (A).

### Option C — Anneal τ_m² (tempering)

Start with inflated `τ_m²` values (so the score `τ_m²·(A_s − B_s)` is
larger by orders of magnitude), iterate to convergence, then gradually
shrink τ_m² to the true prior values while warm-starting. This is a
classical trick for variational inference with weak signals.

Less principled than (A) but easy to implement; works around the
fundamental scaling issue rather than fixing it.

### Recommended path

Start with **(A) Wakefield-style Bayes factor** because:
- It directly addresses the documented failure mode.
- Mirrors the SNP-side update structure — consistent treatment of fixed
  and random effects.
- Has a clean variance-component-LRT interpretation that lines up with
  the manuscript's existing REML framing.
- If (A) is numerically stable at our n, no need for (B) or (C).

Fall back to (B) + (C) if (A) is still too weak.

## Implementation plan (sub-blocks)

### 2.6.5.1 — Derive the Bayes-factor segment update

- [ ] Work out the closed-form (or carefully approximated) Bayes factor
      for `c_s = m` under a single variance-component τ_m² G_s,
      analogous to Wakefield for SNPs.
- [ ] Show its `O(n)` score scaling and document the derivation
      alongside the existing variational derivation in `main.tex`.
- [ ] Smoke-test on a 2-component (spike + single slab) model on a
      tiny sim. Verify the ρ values separate (close to 0 or 1) instead
      of staying near the prior.

### 2.6.5.2 — C implementation

- [ ] New file `src/lib/var_random_mixture.c` (or branch in `var_random.c`)
      with the mixture-prior update path. Existing point-estimate REML EM
      remains as default behaviour.
- [ ] CLI flag `--ws2-prior reml | mixture` (default `reml`).
- [ ] Match the existing data plumbing (Hutchinson trace estimator
      for B_s, CG-based V-solve for r̃).
- [ ] Output `ρ_{s,m}` per tree alongside the legacy `ws²_s = Σ_m ρ_sm τ_m²`.

### 2.6.5.3 — Polygenicity recovery test (the primary target)

Re-run the Block 2.5.2 sweep (3 block_frac × 2 demog × 5 seeds, both
sparse and scattered architectures — per project rule) with
`--ws2-prior mixture`.

- [ ] Compute estimated M_e per cell; compare to truth.
- [ ] **Acceptance #1**: estimated M_e tracks truth M_e within **2×**
      across all 6 cells × both architectures.
- [ ] **Acceptance #2**: variance-capture curves visibly differ at
      bf=0.005 vs bf=0.5 on the estimator side, under both architectures.

### 2.6.5.4 — Prediction-no-regression test (the safety target)

**Per Adam 2026-06-01**: spike-and-slab can hurt prediction by
over-sparsifying. Run the full prediction pipeline (BayesR + ARG-BLUP-uni
+ SARGENT-mixture) on the rare-variant ablation sims:

- [ ] Re-run `rare_variant_ablation.sh` with the new estimator on the
      4 coverages × 2 architectures × 5 seeds = 40 cells.
- [ ] **Acceptance #3**: SARGENT-mixture R² is within ±0.02 of
      SARGENT-REML (current default) on every cell — **no
      systematic prediction regression**. The R² should still track
      ARG-BLUP-uniform.
- [ ] If a regime shows >0.02 R² loss: report it; either tune `τ_m²`
      values to soften the spike or expose a CLI knob (`--mixture-shrinkage`).

### 2.6.5.5 — Localization-no-regression test

The Block 2.6 fused-lasso experiment hurt localization AUC at high λ
(collapsed to chance). Need to verify spike-and-slab doesn't do the
same.

- [ ] Re-run the Block 1.2 localization sweep (30 cells: 3 contrast ×
      2 demog × 5 seeds) with `--ws2-prior mixture`.
- [ ] **Acceptance #4**: estimator AUC on every cell ≥ current (REML)
      estimator's AUC − 0.05 (small tolerance for stochasticity).
- [ ] Ideally AUC *improves* in the cells where REML is currently
      under-performing — the same target the fused-lasso was meant
      to hit.

### 2.6.5.6 — Hyperparameter tuning

- [ ] Choose `τ_m²` values (default mixture vs adaptive).
- [ ] Choose `φ_m` priors (uniform vs sparse).
- [ ] Cross-validation or M-step for `φ_m` (Option B above, if
      Option A alone is insufficient).

### 2.6.5.7 — Write-up

- [ ] `test/eval/block2_6_5_spike_slab.md` with the verdict — analogous
      to `block2_6_verdict.md`.
- [ ] Document the per-segment Bayes-factor derivation in
      `main.tex` Methods, alongside the existing variational derivation
      (which is preserved for completeness — the new path is one
      additional algorithmic option, not a replacement).

## Pre-registered acceptance criteria (consolidated)

| # | criterion | metric | threshold |
|---|---|---|---|
| 1 | Polygenicity recovery | Estimated M_e vs truth M_e on Block 2.5.2 sweep | Within 2× across all 6 cells × 2 archs |
| 2 | Variance-capture distinguishability | Curve shape at bf=0.005 vs bf=0.5 on estimator side | Visible separation, both archs |
| **3** | **No prediction regression** | **SARGENT-mixture R² vs SARGENT-REML R² on rare-variant ablation** | **Within ±0.02 R² on every cell** |
| 4 | No localization regression | Block 1.2 AUC | ≥ REML AUC − 0.05 on every cell |
| 5 | Stable hyperparameter choices | Cross-validated (τ², φ) stable across seeds | CV bands narrow |

**Net verdict requirement**: criteria 1 + 3 must both pass. Criterion 1
alone would help localization but threaten prediction (Adam's caveat);
criterion 3 alone is fine but doesn't justify the implementation effort.

## Effort estimate

- 2.6.5.1 (derivation): 1–2 days. Mostly literature work + math; analogous
  to existing SNP derivation.
- 2.6.5.2 (implementation): 4–7 days. C work in `var_random.c`; parallel
  path to the existing EM, gated by a flag.
- 2.6.5.3 + 2.6.5.4 + 2.6.5.5 (experiments): 1–2 days of compute + analysis.
- 2.6.5.6 (tuning): 2–3 days iterating.
- 2.6.5.7 (writeup): 1 day.

**Total: ~2 weeks.**

## What's at stake

If Block 2.6.5 succeeds (acceptance #1 + #3): the polygenicity-recovery
story (Block 2.5 / Fig 5 in manuscript) is fixed, the variance-capture
curve becomes a manuscript-quality readout, and the spike-and-slab
description in `main.tex` Methods becomes a faithful implementation
description rather than a derivation that the code abandoned.

If 2.6.5 fails acceptance #3 (prediction regresses): we ship it as
opt-in only — same fate as Block 2.6 fused-lasso. Prediction-optimal
estimator stays as the current REML EM (≡ uniform-equivalent for
prediction); spike-and-slab is for localization studies only.

If 2.6.5 fails acceptance #1 (polygenicity still flat): the documented
S>>n scale issue is fundamental; need n in the thousands before the
spike-and-slab can discriminate. The manuscript's polygenicity claim
gets a sample-size caveat.

## Branch + workflow

- Feature branch: `lap/block2.6.5-spike-slab` (will be created when
  implementation starts).
- Default behaviour on `main` is unchanged — `--ws2-prior reml` is the
  default. New behaviour is opt-in via `--ws2-prior mixture`.
- PR against `main` once acceptance criteria pass.

## Open question for Adam before implementation starts

The Wakefield-style Bayes factor (Option A) is the recommended starting
point. **But the derivation needs care** — the SNP analogue treats β_j
as a single scalar; the segment analogue treats ws²_s as a single
scalar variance component but the data are `y` (n-dim) projected onto
`G_s` (rank up to n). The Bayes factor's variance ratio is
`tr(V⁻¹ G_s) · τ_m² + 1` (?), which scales as O(n/S) per tree.

If `n/S · τ_m²` is large enough (e.g., n=300, S=1000, τ_m²=1 → 0.3), the
score becomes meaningful. If still too weak, fall to Option B/C.

**Worth a quick analytic check before committing to (A)**: at typical n
and S, what's `tr(V⁻¹ G_s)`'s magnitude under both the null and a
plausible-active w_s²? If `n/S · τ_m² < 1` in the null vicinity, the BF
won't discriminate either, and we need to think harder. ~30 min of
analytic work + a tiny Python prototype before writing C.
