# Cross-ancestry PRS failure has an allele-age signature (2026-06-12)

**Headline result (F1, the cross-ancestry flagship — now decisive):**

> The portable component of heritability across ancestries is the **deep,
> shared-ancestral (old) variation**; the recent, population-specific (young)
> variation does not transfer. The ARG's **allele-age** axis identifies the
> portable component, and **frequency cannot** — once you know a variant's age,
> its frequency tells you nothing more about whether it transfers. And
> cross-ancestry prediction is **capped by the ancient/shared heritability
> fraction**, independent of within-ancestry power — a mechanistic explanation
> for real-world cross-ancestry PRS collapse.

This is **prediction**, not variance decomposition — so it sidesteps the
collinearity wall that killed the per-bin GREML partitioning (see
`block2_2_design.md` Phase C verdict + ROADMAP 2026-06-09). Age is genuinely
ARG-distinctive (genotypes can't date variants), and the result is tied to
SARGENT's one confirmed strength (Block 1.1 ARG-BLUP cross-ancestry portability).

Setup throughout: `--demography transfer`, pure A→B (train on pop A, test on
pop B, `mix_fraction_B=0`), two split depths F_ST≈0.15 (split 3500) and ≈0.30
(7000). "old/young" = age (mutation time) median split; "common/rare" =
MAF≥/<0.05; eGRM age bins via `kinship_assembly.bin_by_age_tertile`.

## 1. Mechanism — cross-ancestry relatedness is 100% in the deep stratum

The train-A × test-B cross-kinship, decomposed by age stratum (no phenotype, no
fitting — pure genealogy), F_ST≈0.30:

| stratum | cross-kinship (A×B) | % of baseline | cross/within |
|---|---|---|---|
| recent | 0 | 0% | 0.00 |
| middle | 0 | 0% | 0.00 |
| **ancient** | **= baseline** | **100%** | **0.66** |

A and B are related **only** through pre-split branches. Cross-ancestry
prediction can ride **only** on this shared-ancestral relatedness.

## 2. Decisive test — age predicts transfer, frequency does not (z=15–20)

Per-variant transfer prerequisite `T = var_B/var_A` (does the variant carry
variance in the target pop), regressed on age + frequency; partial correlations
meta-analyzed across 30 seeds × 2 depths (~100k variants):

| | F_ST≈0.15 | F_ST≈0.30 |
|---|---|---|
| **partial corr(AGE \| freq)** | **0.139, z=19.7** | **0.107, z=15.0** |
| partial corr(FREQ \| age) | −0.020, z=−2.0 | 0.011, z=1.5 |
| ΔR² adding age to a freq-only model | 0.021 | 0.013 |

Age predicts transferability strongly controlling for frequency; **frequency,
controlling for age, adds essentially nothing.** r²-independent (genotype level),
so immune to the prediction-R² issue below.

Visualized — fraction of A-variants also present in B, by age × MAF tertile
(pooled, F_ST≈0.30):

| | MAF rare | MAF mid | MAF common |
|---|---|---|---|
| **AGE ancient** | 42% | 40% | 64% |
| AGE middle | 0.8% | 0.3% | 0.4% |
| AGE recent | 0.04% | 0.07% | 0% |

Sharing drops ~3 orders of magnitude **down the age rows** and is ~flat
**across the frequency columns** — the heatmap stripes horizontally (age), not
vertically (frequency).

## 3. Prediction-level — portability stratifies by age (and frequency)

ARG-BLUP portability = cross-R²(A→B) / within-R²(A→A), 30 seeds:

| contrast | F_ST≈0.15 | F_ST≈0.30 |
|---|---|---|
| FREQ (common > rare) | +0.119, z=2.3 | +0.045, z=2.7 |
| **AGE (old > young)** | **+0.107, z=3.0** | **+0.052, z=3.2** |
| BEYOND (old_rare > young_rare) | +0.054, z=1.4 | +0.020, z=1.9 |

Old/shared ports ~4× young/specific; the ARG **reproduces the frequency story**
(common>rare) **and adds age on top** (age is the cleaner stratifier, z=3 vs
2.5). The strict within-frequency age residual (old_rare vs young_rare) is
suggestive but not conclusive at the prediction level (z=1.4–1.9) — it is the
genotype-level regression (z=15–20) that nails "age beyond frequency".

## 4. The cap — portability is bounded by the shared fraction

Re-running the BLUP at higher power (n=2800, h²=0.8, less polygenic): within-R²
rose to ~0.5, but **cross-R² stayed pinned at ~0.01–0.02** — portability ~4% at
F_ST≈0.30, the same as at lower power.

> Cross-ancestry prediction does **not** scale with within-ancestry power. It is
> capped by the small portable (ancient/shared) heritability fraction. You can
> predict within-A arbitrarily well; A→B is bounded by what is shared.

This is the mechanistic statement of real-world cross-ancestry PRS collapse, and
it means the BLUP *portability ratio* cannot be "powered up" by raising r² — the
genotype-level regression (§2) is the high-power test, not the ratio.

## Proposed figure — "Cross-ancestry PRS failure has an allele-age signature"

- **Panel A (mechanism/cartoon):** two populations splitting; deep ancestral
  branches shared (carry prediction), recent branches population-specific (don't).
  Inset bar: cross-kinship by age stratum (recent 0 / middle 0 / ancient 100%).
- **Panel B (the decisive visual):** the age × MAF heatmap of "fraction shared
  with B" — horizontal stripes (age), flat across frequency. Annotate partial
  corr(age|freq) z=15–20 vs (freq|age) z≈0. *This is the money panel.*
- **Panel C (prediction payoff):** portability (cross/within R²) by age stratum
  (old ports, young doesn't), with the frequency split alongside — "ARG recovers
  the frequency baseline and adds age."
- **Panel D (the cap):** within-R² vs cross-R² across regimes — within rises to
  0.5, cross flat at 0.02. Caption: portability bounded by the shared fraction.

**Integration with current figures.** This is the **mechanistic companion to
the existing cross-ancestry portability figure** (Block 1.1: ARG-BLUP > typed-SNP
cross-ancestry) — it *explains why* ARG-BLUP ports better (it captures the
shared/ancient component). It is the natural **F1 flagship figure**, and it can
**fill the slot vacated by the dropped Fig 5 Panel A** (functional "WHERE",
retired per the 2026-06-09 identifiability audit). Suggest placing it
immediately after the Block 1.1 portability figure as its "why".

## Reproduce

```bash
PYBIN=/DATA/home/asiepel/miniforge3/envs/sargent/bin/python
# sims: transfer A->B, e.g. F_ST~0.30:
$PYBIN scripts/simulate.py --demography transfer --transfer-split-time 7000 \
  --mix-fraction-B 0.0 --n-indiv 1000 --n-test 400 --length 1e6 --Ne 1e4 \
  --h2-snp 0.5 --h2-arg 0.0 --causal-fraction 0.02 --seed 1 --out <cell>/case
# decisive regression (age vs freq for transfer):  test/eval/xanc_transfer_regress.py
# BLUP portability by stratum:                      test/eval/xanc_portability.py
#   (CELL=<cell> env var; reads case.full.{trees,phen,params.json})
```
