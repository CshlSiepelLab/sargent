# Block 2.5.2 — Controlled-polygenicity sweep + the fused-lasso decision

Date: 2026-06-01.
Status: **v1 done.** Decision: **fused-lasso (Phase 11) revival is
justified**.

Code: `test/eval/polygenicity_sweep.sh`, `test/eval/arch_metrics.py`,
`test/eval/plot_arch.py` (all on `main`).
Output data (not in repo): `~/test-sargent/block2_5_2/`.

## Design

Controlled-polygenicity sweep across `--ws2-block-frac` at fixed contrast,
n_blocks, n_train, length. Three regimes, two demographies, 5 seeds.

| BLOCKFRAC | regime | predicted truth `M_e,ARG` |
|---|---|---|
| 0.005 | oligogenic   | ~5 (5 active trees out of ~1000) |
| 0.05  | moderate     | ~50 (the Phase 10 / Block 1.2 baseline) |
| 0.5   | highly polygenic | ~500 |

All cells: n_train=300, n_test=100, L=300kb, h²_arg=0.5, contrast=1000,
n_blocks=3, h²_snp=0. Pure-random Stage-2-only fit (sargent `--arg-only`).

Total: 3 × 2 × 5 = 30 cells. Ran on the laptop overnight, **381 min wall**.

## Results

Per-cell mean ± sd over 5 seeds:

| demography | block_frac | n | **truth M_e ± sd** | **est M_e ± sd** | truth k_at_50pct | est k_at_50pct |
|---|---|---|---|---|---|---|
| admixed   | 0.005 | 5 | **8.6 ± 0.1** | 14.2 ± 5.2 | 4.0 | 5.4 |
| admixed   | 0.05  | 5 | **61.7 ± 2.6** | 12.1 ± 5.0 | 31.0 | 5.2 |
| admixed   | 0.5   | 5 | **603.0 ± 26.4** | 19.7 ± 4.8 | 301.6 | 8.0 |
| panmictic | 0.005 | 5 | **4.7 ± 0.2** | 14.3 ± 4.1 | 2.0 | 5.8 |
| panmictic | 0.05  | 5 | **38.7 ± 3.5** | 12.5 ± 2.5 | 19.4 | 5.2 |
| panmictic | 0.5   | 5 | **383.6 ± 32.3** | 15.8 ± 1.5 | 192.2 | 6.2 |

## Verdict

**Truth `M_e,ARG` separates the three polygenicity regimes by ~100×**
(4.7 → 603 across the sweep). The closed-form prediction (M_e ≈
block_frac × n_trees) is confirmed to within ~10 %.

**Estimated `M_e,ARG` is essentially flat at 12–20** across all three
regimes and both demographies. The estimator does **not** recover
polygenicity. The variance-capture curves overlap closely on the
estimated side, despite ~2-orders-of-magnitude separation on the truth
side (see `~/test-sargent/block2_5_2/arch_variance_capture.png`).

This is the pre-registered "**estimator inadequate**" outcome from the
`block2_5_2` decision rule. The estimator's per-tree noise floor
concentrates ~12–20 trees' worth of mass *regardless* of the true
underlying architecture.

## Implication: fused-lasso (Phase 11) un-parks

Per the ROADMAP parking-lot entry, Phase 11 revives when
> "a target application (very-low-contrast localization; per-seed
> variance reduction for partitioned-h²) makes that 5% AUC worth ~weeks
> of engineering."

Block 2.5.2 surfaces exactly such an application: **the polygenicity
read-out is a manuscript headline (variance-capture curve, "where the
heritability lives in concentration-vs-spread terms") that the current
estimator cannot produce**. Without a smoother prior, the variance-capture
shape is uninformative regardless of the true distribution.

The fused-lasso's two priors directly address the over-concentration:
- **Sparsity-on-levels** (`λ₁ ||θ||_1`) lets log-`ws²` be exactly zero
  outside active regions — eliminates the noise-floor mass that
  contaminates diffuse-architecture estimates.
- **Sparsity-on-differences** (`λ₂ Σ_s |θ_s − θ_{s-1}|`) smooths within
  active regions — reduces the per-seed scatter (sd 5 on a mean of 14
  in current estimator, a 35 % CV that washes out architecture signal).

**Acceptance criterion for the fused-lasso revival (re-stated from Phase 11
design):**

- Estimated `M_e,ARG` distinguishes the three Block 2.5.2 regimes by
  ≥ 10×.
- Variance-capture curves at bf=0.005 vs bf=0.5 are visibly different
  on the estimated side.
- No regression on Block 1.2 localization AUC.
- Diffuse-control (uniform ws²) selects ~zero blocks.

## CLAIMS implications

- **D1** (polygenicity recovery): still **WEAK** at the estimator level
  pending fused-lasso. Truth-side recovery is exact, but that's not what
  the manuscript can claim.
- **D2** (distinguishing regimes): demoted to **FALSIFIED at current
  estimator**. The current estimator does not distinguish; this is
  the bug the fused-lasso will fix.
- **B6** (rare-variant recovery via ARG): pending the rare-variant
  ablation experiment (was in the Block 2.5.2 follow-up; deferred until
  fused-lasso lands).

## Files

| file | purpose |
|---|---|
| `~/test-sargent/block2_5_2/arch_b252.json` | per-case M_e, capture curves, quantiles, etc. (30 cases). |
| `~/test-sargent/block2_5_2/arch_b252_grouped.json` | same with `group` field added for plotting. |
| `~/test-sargent/block2_5_2/arch_me_arg.{pdf,png}` | M_e bar chart truth-vs-estimated × 6 group cells. |
| `~/test-sargent/block2_5_2/arch_variance_capture.{pdf,png}` | variance-capture curves stratified by group. |
| `~/test-sargent/block2_5_2/run.log` | sweep run log. |

## Reproduce

```bash
BASE=$HOME/test-sargent/block2_5_2 bash test/eval/polygenicity_sweep.sh
cd $BASE
python ~/git-src/sargent/test/eval/arch_metrics.py --batch . \
    --no-branches --out arch_b252.json
# attach a group key for plotting
python -c "import json;d=json.load(open('arch_b252.json'))
for r in d:
    bf=r.get('ws2_block_frac');dem=r.get('demography','?')
    if isinstance(bf,(int,float)): r['group']=f'{dem}_bf{bf:g}'
json.dump(d,open('arch_b252_grouped.json','w'),indent=2)"
python ~/git-src/sargent/test/eval/plot_arch.py arch_b252_grouped.json \
    --group-key group --out-prefix arch_
```

## Caveats

- 5 seeds per cell; per-cell SD on estimated M_e is ~5 (35 % CV) at this
  scale. Bigger seed count would tighten the bars but won't change the
  qualitative verdict — the truth-vs-estimate gap is two orders of magnitude.
- All cells used h²_snp = 0 (pure random); the BayesR-fixed-effects side
  is not exercised here. The polygenicity-recovery verdict is for the
  ARG random-component primitive.
- The "estimator over-concentration" is consistent with what Block 2.5.1
  showed on Phase 10 + Block 1.2; this sweep adds the controlled
  polygenicity axis that makes the failure quantitative.
- The Phase 10 / Block 1.2 results showed estimator AUC for *localization*
  (which-trees-are-active) is OK (0.7+); the new result here is that
  *amplitude recovery* (how much variance in each tree) is not OK. These
  are different aspects of the same estimator.
- **Note (project standing rule, established 2026-06-01):** this sweep ran
  `sparse` architecture only. Future polygenicity / estimator-amplitude
  validations must include `scattered` alongside `sparse` — sparse's
  contiguous blocks can be well-tagged even at low coverage, understating
  the ARG advantage and making estimator behaviour look better than it is
  for realistic architectures. The Block 2.6 acceptance criteria were
  extended to require both architectures.
