# Phase 10 — Oracle diagnostic: can per-segment `ws²` ever beat ARG-BLUP-uniform?

Date: 2026-05-28.
Code: `scripts/simulate.py` (`--ws2-architecture`, `--demography`),
`test/eval/oracle_ws2.py`, `test/eval/oracle_ws2.sh`.
Data: `~/test-sargent/phase10_admix/`.

## Design

A purely-random simulation (`--h2-snp 0`, `--h2-arg 0.5`) under genuinely
favourable conditions for per-segment weights, comparing three paired
ARG-BLUP arms that differ only in their per-tree weight vector:

- **uniform**: `wₛ = 1/S` (ARG-BLUP-uniform).
- **oracle**: `wₛ = true wsₛ²` from the simulation (`.full.true_ws2.tsv`).
- **estimated**: `wₛ = sargent ŵ²` from `--arg-only` (`.w2.tsv`).

**Setting:**
- **Admixed demography**, two source pops of `Nₑ=10⁴` split 6000 generations
  ago (F_ST ≈ 0.26), recent admixture pulse at 20 generations, α=0.5. Local
  ancestry varies sharply along the genome → different regions carry
  genuinely different relatedness, the regime where weight heterogeneity
  can matter at all.
- **Sparse architecture**: 3 contiguous large-effect blocks spanning 5% of
  trees, block:floor weight ratio 1000 (nearly all variance concentrated in
  the active blocks). Plus a **diffuse control** (flat constant `wₛ²`).
- n_train = 300, n_test = 100, L = 300 kb (~700–1000 trees), 5 seeds per arm.

**Two ceilings measured:**

1. **Prediction ceiling** = `R²_oracle − R²_uniform`. Upper bound on what
   *any* per-segment estimator could deliver in OOS prediction over uniform.
2. **Localization ceiling** = AUC of the idealized per-tree score
   `Aₛ / tr(Gₛ)` evaluated at the homogeneous-null `V_uniform` (exact, no
   Hutchinson/estimation noise; computed in `idealized_loc_auc`). Upper
   bound on any score-based localizer's active-block detection.

## Results

Full per-seed table (also in `phase10_admix/phase10_table.md`):

| arch | seed | R²_uni | R²_oracle | R²_est | **oracle−uni** | est−uni | **loc-ceil AUC** | est AUC |
|---|---|---|---|---|---|---|---|---|
| diffuse | 1 | +0.458 | +0.458 | +0.445 | +0.000 | −0.013 | — | — |
| diffuse | 2 | +0.200 | +0.200 | +0.125 | +0.000 | −0.076 | — | — |
| diffuse | 3 | +0.532 | +0.532 | +0.514 | +0.000 | −0.018 | — | — |
| diffuse | 4 | +0.418 | +0.418 | +0.412 | +0.000 | −0.006 | — | — |
| diffuse | 5 | +0.384 | +0.384 | +0.367 | +0.000 | −0.017 | — | — |
| sparse  | 1 | +0.322 | +0.356 | +0.335 | **+0.034** | +0.013 | **0.843** | 0.902 |
| sparse  | 2 | +0.332 | +0.361 | +0.325 | **+0.029** | −0.006 | **0.904** | 0.729 |
| sparse  | 3 | +0.464 | +0.484 | +0.442 | **+0.020** | −0.022 | **0.778** | 0.594 |
| sparse  | 4 | +0.376 | +0.407 | +0.383 | **+0.031** | +0.007 | **0.856** | 0.732 |
| sparse  | 5 | +0.388 | +0.389 | +0.356 | **+0.001** | −0.032 | **0.837** | 0.743 |

Summary over the 5 sparse seeds (h²_arg realised 0.506, target 0.5):

| metric | mean | sd | min | max |
|---|---|---|---|---|
| Prediction ceiling (oracle − uniform R²) | **+0.023** | 0.013 | +0.001 | +0.034 |
| Estimator gap (estimated − uniform R²)   | −0.008 | 0.019 | −0.032 | +0.013 |
| Localization ceiling (idealized AUC)     | **0.844** | 0.045 | 0.778 | 0.904 |
| Estimator localization AUC               | 0.740 | 0.109 | 0.594 | 0.902 |

Paired SE on the prediction-ceiling Δ ≈ 0.006 (n=5); 95% CI ≈ +0.01 to +0.04.
All five sparse seeds positive and tightly clustered. Diffuse-control
oracle = uniform to the digit on all 5 seeds (pipeline sanity ✓).

## Interpretation

The two ceilings answer the two questions in opposite directions:

**Prediction ceiling is firm and small.** Even under deep admixture
(F_ST≈0.26) with concentrated 5%-active blocks, knowing the *true* per-segment
weights buys only **+0.023 R²** over ARG-BLUP-uniform. The current per-tree
REML captures none of it (slightly worse than uniform on average). So
per-segment `ws²` heterogeneity is **not** a meaningful prediction lever
even with a perfect estimator — the ARG kinship + uniform weighting already
captures the bulk of the prediction signal. The earlier panmictic pilot
gave +0.046 from one seed; the rigorous admixed 5-seed mean is *smaller*,
not larger, confirming that the prediction ceiling is structurally bounded
near +0.02 R² and that adding favourable structure does not lift it.

**Localization potential is substantial.** Idealized AUC averages 0.84
under structure — active blocks are clearly identifiable in principle. The
current estimator reaches AUC 0.74 mean (already well above chance, and
much better than the AUC ≈ 0.5 from panmictic small-scale debugging), but
with large per-seed variance (0.59–0.90) and a robust ~0.10 gap below the
ceiling. **There is real headroom for a better estimator** — closing
0.74 → 0.84 in localization AUC.

## Verdict

> **Per-segment `ws²` is a localization tool, not a prediction tool.**

This separates Phase 11's scope cleanly. The fused-lasso redesign should be
pursued as a **localization / local-`h²` method**, not as a prediction
enhancement. Its self-certifying degradation to ARG-BLUP-uniform (`λ→∞`)
protects prediction (it can't lose materially given the bounded ceiling)
while the sparsity-on-levels + sparsity-on-differences prior, together with
ARG-aware fusion weights, directly attack the localization gap that has
room to move (0.74 → 0.84+ AUC).

The manuscript's *prediction* story belongs to the ARG kinship + the right
demography (transfer / admixture / Neandertal), where ARG-BLUP-uniform is
the workhorse. The *localization* / partitioned-`h²` story is where the
per-segment-weights contribution earns its keep.

## Caveats and what was not tested

- 5 seeds is modest. CIs are wide but the all-positive, tightly-clustered
  pattern is unambiguous on the qualitative verdict.
- Even more extreme architectures (1% active, contrast 10⁴, single mega-locus)
  might lift the prediction ceiling further, but at 5%/1000 we are already in
  an unrealistically favourable regime; the bound is unlikely to grow much.
- Active blocks were placed at evenly-spaced positions, not aligned with
  the genome's most-structured regions. Aligning blocks with maximum
  ancestry-divergent windows could push both ceilings up modestly.
- The "localization ceiling" uses the heuristic `Aₛ/tr(Gₛ)` detector at
  exact `V_uniform`; the proper REML score `Aₛ − Bₛ` is slightly tighter
  but requires `tr(V⁻¹Gₛ)` per branch (expensive). In 1/5 seeds the
  estimator exceeded this ceiling (s1), which just means the heuristic
  detector is not strictly optimal; the qualitative ordering
  (estimator ≲ ceiling, ceiling ≫ chance) holds.

## Reproduce

```bash
DEMOG=admixed SPLIT=6000 ADMIXT=20 ALPHA=0.5 \
  N_TRAIN=300 N_TEST=100 LENGTH=300000 H2_ARG=0.5 \
  BLOCKFRAC=0.05 NBLOCKS=3 CONTRAST=1000 \
  SEEDS="1 2 3 4 5" ARCHS="diffuse sparse" \
  BASE=/path/to/phase10_admix \
  bash test/eval/oracle_ws2.sh
```

Outputs land per-case at `<BASE>/<arch>_s<seed>/est.oracle.json`; the driver
also writes `<BASE>/phase10_table.md`.
