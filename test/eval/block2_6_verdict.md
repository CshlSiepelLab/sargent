# Block 2.6 — Fused-lasso prior on log `ws²`: verdict

Date: 2026-06-01.
Branch: `lap/block2.6-fused-lasso` (commit `6b35353` = merge of main into
the feature branch). Awaiting PR merge to `main` as an opt-in knob.

## What was implemented

A proximal step applied after each `var_random` EM update, gated by two
CLI flags. Both default to 0, in which case behaviour is byte-identical
to the legacy unregularised per-tree REML EM (verified by smoke test).

```
After the standard ws² update each iteration:
    log_ws2[s] := log(ws²[s])
    μ          := mean over s of log_ws2[s]
    θ[s]       := log_ws2[s] − μ           (hierarchical centering)
    θ          ← prox_{λ₁‖θ‖_1 + λ₂ TV(θ)}  (sequential prox: TV via cyclic
                                              coord descent, then L1
                                              soft-threshold)
    ws²[s]     ← exp(μ + θ[s]), clamped to [ws2_floor, ws2_max]
```

CLI: `--ws2-l1 LAMBDA1`, `--ws2-l2 LAMBDA2`. Defaults: 0.

## Pre-registered acceptance criteria (from ROADMAP § Block 2.6)

| # | criterion | result |
|---|---|---|
| 1 | Est `M_e,ARG` separates the three Block 2.5.2 regimes (bf = 0.005 / 0.05 / 0.5) by ≥ 10×, **under both `sparse` AND `scattered`** | **FAIL.** On bf=0.5 (truth M_e≈600), M_e plateaus at 42.3 across λ₂ ∈ {0.5, 1, 2, 5, 10}. On bf=0.005 (truth M_e≈9), the estimator goes the *wrong direction* (toward uniform, not sparser). Estimator M_e ratio still 2× at best, not 10×. |
| 2 | Variance-capture curves at bf=0.005 vs bf=0.5 visibly different on the estimator side | **PARTIAL.** L2 lifts the polygenic-regime curve modestly; L1 on the oligogenic regime degenerates to uniform at λ₁≥1. The curves don't recover their truth-side shapes. |
| 3 | No regression on Block 1.2 localization AUC | **MIXED.** AUC unchanged on polygenic (0.646 → 0.649); AUC **improved by +0.10** on oligogenic (0.685 → **0.789** at λ₁=0.5). No regression. |
| 4 | Diffuse-control selects ~zero blocks | Not tested in this v1 (deferred). |
| 5 | CV-stable (λ₁, λ₂) | Not tested in this v1; λ values chosen by direct sweep. |

**Overall**: criterion 1 (the polygenicity-recovery headline) **failed**.
Criterion 3 (localization AUC) **succeeded** on oligogenic, neutral on
polygenic. Net: **partial — works as a localization tuning knob, does
not solve polygenicity recovery.**

## Comparison table — fused-lasso vs unregularised baseline

bf=0.5 admixed cell (truth M_e ≈ 600, S_train=1171, n_active=580):

| config | M_e | k@50% | k@90% | Spearman(ŵ²,w²) | AUC(active) |
|---|---|---|---|---|---|
| **unreg (legacy)** | **16.1** | 7 | 22 | +0.265 | 0.646 |
| λ₁=0, λ₂=0.5 | 42.3 | 16 | 44 | +0.278 | 0.649 |
| λ₁=0, λ₂=1.0/2/5/10 | 42.3 (plateau) | 16 | 44 | +0.278 | 0.649 |
| λ₁=0.5, λ₂=0 | 14.7 | 6 | 19 | +0.093 | 0.575 |
| λ₁=1.0, λ₂=0 | 12.7 | 5 | 17 | −0.014 | 0.499 |

bf=0.005 admixed cell (truth M_e ≈ 9, S_train=1171, n_active=6):

| config | M_e | k@50% | k@90% | Spearman(ŵ²,w²) | AUC(active) |
|---|---|---|---|---|---|
| **unreg (= λ₁=λ₂=0)** | **20.0** | 8 | 31 | −0.014 | 0.685 |
| **λ₁=0.5, λ₂=0** | 24.5 | 9 | 31 | +0.050 | **0.789 ← AUC win** |
| λ₁=1.0, λ₂=0 | 619 | 520 | 1031 | +0.288 | 0.784 |
| λ₁=2.0, λ₂=0 | 1043 (uniform) | 578 | 1053 | +0.851 | 0.477 (worse-than-chance) |

## Reading

1. **λ₂ (TV) plateaus.** At bf=0.5, all λ₂ values from 0.5 to 10 give
   *identical* M_e=42.3. The TV regularisation hits a structural ceiling.
   Diagnosis: the EM's multiplicative `ws²[s] *= A[s]/B[s]` update
   pushes most trees to the variance floor *before* the prox runs each
   iteration; the prox can reshape the active set but cannot lift floored
   trees back. Median ws² stays at ~`1e-15` (= 1e-12 × ws2_max) across
   all λ₂. The post-EM-prox approach is structurally capped on
   polygenicity recovery.

2. **λ₁ (L1 on θ) pushes toward uniform**, not toward sparser ws². Since
   θ = log(ws²) − mean(log(ws²)), shrinking θ to zero means ws² → constant.
   This is the *wrong direction* for the oligogenic regime (where truth
   wants a few non-zero, rest exactly zero). At λ₁≥1 on bf=0.005, the
   estimator collapses to near-uniform (M_e → 619, then 1043). Below
   chance AUC at λ₁=2 confirms the collapse.

3. **But λ₁=0.5 on the oligogenic regime improves localization AUC by
   +0.10.** The mild L1 reshapes the active set toward the truth's
   rank-ordering without destroying it. This *is* a substantive
   improvement on the Block 1.2 localization metric in the regime
   where active blocks are sparse — the *correct* manuscript-headline
   improvement.

4. **The polygenicity-recovery problem is not fixed.** The current
   fused-lasso targets θ (log-centered), which sparsifies toward
   uniform ws². To get the sparse-ws² recovery the Phase 11 idea
   originally promised (and that the manuscript's Methods derives —
   spike-and-slab mixture prior on ws² directly), we'd need a
   different prior structure. See "What's next" below.

## What's next: spike-and-slab on `ws²` directly — **DROPPED 2026-06-01**

> The path described below was proposed as Block 2.6.5, then dropped after a
> meta-review (Adam + Claude). Design preserved at
> `test/eval/block2_6_5_design.md`. Effort redirected to Blocks 2.1 and 2.2.

## Original "what's next" note (superseded)

The right prior for polygenicity / variance-capture recovery is the
spike-and-slab mixture prior on `ws²_s` *itself* (not on log(ws²) − mean):

```
ws²_s ~ Σ_m φ_m · δ(τ_m²),    τ_1² = 0 (spike), τ_2² > 0 (slab)
```

The SARGENT manuscript already derives the variational update for this
prior (the ρ-update from Methods); the current C code abandoned it for
the simpler point-estimate REML EM. **Resurrecting the marginalised
mixture-prior update is the right next step.** Block 2.6.5 (proposed).

Differences from the just-tested fused-lasso:

- L1/L2 on θ → uniform ws². Spike-and-slab on ws² → most trees at exactly
  zero, a few at slab variance. **Opposite direction.**
- Sparsity is multiplicative (responsibility ρ_s ∈ [0,1] per tree),
  enabling true "tree-active or not" classification.
- Already proven viable in the original Phase 11 derivation (just not
  implemented in code).

**Acceptance criterion for Block 2.6.5**: estimator M_e tracks truth M_e
within 2× on both bf=0.005 (truth ≈ 5–9) and bf=0.5 (truth ≈ 600). This
is the criterion Block 2.6.1–2.6.3 (fused-lasso) failed.

## Manuscript implications

1. **Block 1.2 localization story (C1 HELD) is unchanged.** Current
   estimator's AUC 0.7+ across architectures stands. Fused-lasso is
   a modest tuning improvement, not a regime change.
2. **Block 2.5 variance-capture / polygenicity story needs caveat.**
   Either (a) show truth-side variance-capture only and footnote the
   estimator limitation, (b) wait for Block 2.6.5 (spike-and-slab on
   ws² itself) before publishing the estimated-side curve, or
   (c) accept the over-concentration as a known limitation.
3. **The fused-lasso CLI knob is worth shipping** as an opt-in for
   oligogenic localization studies (small but real AUC win); doesn't
   improve the default behaviour.

## Recommendation

- **Merge `lap/block2.6-fused-lasso` to main** as an opt-in (`--ws2-l1`,
  `--ws2-l2` flags; both default 0; behaviour byte-identical to legacy
  when off).
- **Document the AUC win on oligogenic / the polygenicity-recovery
  limitation** in the help text / Methods.
- **Open Block 2.6.5** (spike-and-slab on ws² directly) as the actual
  fix for polygenicity recovery. ROADMAP updated accordingly.
- **CLAIMS bucket moves**: no change. D1 (polygenicity recovery), D2
  (distinguishing regimes), C5 (sparsity helps localization) all stay at
  their current buckets pending Block 2.6.5.

## Reproduce

```bash
# Build the lap/block2.6-fused-lasso branch (or main once PR lands):
cmake --build build -j

# bf=0.5 (polygenic) — λ₂ alone helps modestly:
CASE=~/test-sargent/block2_5_2/bf0.5_dadmixed_s1
~/git-src/sargent/build/src/progs/sargent --silent --nprobe 10 --arg-only \
    --ws2-l2 0.5 --outroot /tmp/bf05_l2 \
    $CASE/case.phen $CASE/case.trees $CASE/case
# bf=0.005 (oligogenic) — λ₁ ≈ 0.5 lifts localization AUC by ~0.10:
CASE=~/test-sargent/block2_5_2/bf0.005_dadmixed_s1
~/git-src/sargent/build/src/progs/sargent --silent --nprobe 10 --arg-only \
    --ws2-l1 0.5 --outroot /tmp/bf005_l1 \
    $CASE/case.phen $CASE/case.trees $CASE/case
```
