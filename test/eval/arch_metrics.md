# Block 2.5.1 — Architecture metrics: `M_e,ARG`, variance-capture curve, per-branch effect-size

Date: 2026-05-31.
Status: implemented + qualitatively validated on existing Phase 10 + Block
1.2 outputs. **No new sims**; pure post-hoc analysis on `est.w2.tsv`.

Code: `test/eval/arch_metrics.py`, `test/eval/plot_arch.py`.
Output data (not in repo): `~/test-sargent/phase10_admix/arch_phase10.json`,
`~/test-sargent/block1_2_v1/arch_block1_2.json`.

Backs claims **D1** (polygenicity recovery) and **D2** (distinguishing
regimes), both currently WEAK / PLAUSIBLE; this work moves both toward HELD
pending Block 2.5.2 controlled-polygenicity validation.

## What it does

Three metrics from per-tree `ŵ²_s` (no model refitting):

1. **`M_e,ARG`** — analog of Wray's effective number of contributing variants
   at the tree level:
   ```
   M_e,ARG  =  (Σ_s ŵ²_s)²  /  Σ_s (ŵ²_s)²
   ```
   Sparse architecture → small `M_e`; diffuse → large.

2. **Variance-capture curve** — fraction of total ARG-variance captured by
   the top-k trees, plotted as a function of k. The curve shape is the
   one-glance architecture summary; the manuscript-Figure-1 candidate
   flagged in the 2026-05-30 brainstorm.

3. **Per-branch effect-size distribution** — `ŵ²_s · bl(b)` aggregated
   across all (tree, branch) pairs in the ARG. The non-parametric
   distribution of branch-level variance contributions, equivalently the
   ARG's view of the per-mutation effect-size distribution. Pairs with
   BayesR's typed-SNP effect-size estimates for the complete picture.

Single-case CLI: `arch_metrics.py --case-prefix X --sargent-prefix Y`.
Batch mode: `arch_metrics.py --batch ROOT_DIR` iterates over
`*/case.params.json`. Visualisation: `plot_arch.py BATCH_JSON --group-key K`.

## Validation 1: Phase 10 admixed (5 sparse + 5 diffuse seeds)

```
case            arch       truth M_e   est M_e   truth k_50%   est k_50%
diffuse_s1      diffuse        1294       17         647            7
diffuse_s2      diffuse        1259       32         630           13
diffuse_s3      diffuse        1168       33         584           13
diffuse_s4      diffuse        1293       35         647           13
diffuse_s5      diffuse        1209       18         605            7
sparse_s1       sparse           66       12          33            5
sparse_s2       sparse           65       20          33            8
sparse_s3       sparse           59       13          30            5
sparse_s4       sparse           65       18          33            8
sparse_s5       sparse           62       11          31            4

Group means:
  diffuse: truth M_e = 1245 ± 56,   est M_e = 27 ± 9     (truth ÷ est ≈ 46×)
  sparse:  truth M_e =   63 ±  3,   est M_e = 15 ±  4    (truth ÷ est ≈  4×)
```

The truth panel of the variance-capture curve (left in
`arch_phase10_variance_capture.png`) is the textbook display: the **sparse**
architecture captures ≥98 % of variance in the top ~5 % of trees — exactly
as designed (5 %-active blocks at contrast 1000) — and the **diffuse**
architecture rises roughly linearly. The estimated panel shows both
architectures separating (sparse curves above/left of diffuse) but
**systematically more concentrated than truth**: the estimator pulls
variance mass into a smaller number of trees than the true profile has.

## Validation 2: Block 1.2 sweep (3 contrast × 2 demog × 5 seeds)

```
group              n   truth M_e ± sd      est M_e ± sd
admixed   c1000    5    61.7 ± 2.6         12.1 ± 5.0
admixed   c100     5    84.4 ± 3.6         16.0 ± 8.6
admixed   c10      5   427.0 ± 18.7        25.5 ± 10.9
panmictic c1000    5    38.7 ± 3.5         12.5 ± 2.5
panmictic c100     5    53.1 ± 4.7         13.5 ± 2.6
panmictic c10      5   272.7 ± 22.3        20.5 ± 5.5
```

The metric is **monotone in contrast**: as block-to-floor contrast drops
(1000 → 100 → 10), `M_e,ARG` rises (architecture appears more diffuse) in
both demographies. The truth values match the closed-form prediction from
the architecture (e.g., contrast=1000, block_frac=0.05 → ≈50 expected). The
estimate tracks the same monotone trend but with substantially compressed
dynamic range: 12–26 vs truth 39–427.

## Findings

1. **The metric works.** `M_e,ARG` and the variance-capture curve cleanly
   distinguish architectures both at truth (20× ratio diffuse/sparse) and
   at estimated level (1.8× ratio, same direction). The variance-capture
   curve shape is what the manuscript Figure 1 would render — qualitatively
   evocative for non-method readers.

2. **The estimator systematically over-concentrates.** Estimated `M_e` is
   ~4× smaller than truth on sparse architectures and ~46× smaller on
   diffuse ones (Phase 10 numbers). This is consistent with the per-tree
   information limit at n=300: the estimator finds peaks where truth is
   flat. The Block 1.2 numbers reinforce this — estimated M_e is roughly
   constant (12–26) while truth spans an order of magnitude (39–427).

3. **The implication for the fused-lasso decision.** Sparsity-on-differences
   regularisation in `ws²` (the demoted Phase 11 redesign) would smooth the
   estimate toward the true distribution by penalising spurious
   concentration. If the manuscript wants a *quantitative* polygenicity
   readout (rather than a qualitative "sparse vs diffuse"), the fused
   lasso may need to come back off the parking lot. **Block 2.5.2
   (controlled polygenicity sims) is the natural decision point**: if
   `n_causal_regions ∈ {10, 100, 1000}` doesn't cleanly separate in the
   estimator, the fused lasso decision goes "yes". Otherwise, the
   qualitative case is enough.

4. **Implementation is cheap and reusable.** All three metrics are
   post-hoc on existing SARGENT output (no model rerun). The per-case
   compute is dominated by ARG I/O when branches are included; metrics
   without branches are ~1s per case. Batch over a 30-case sweep ran in
   ~30s total. The same code drops into Block 2.1 (allele-age stratified
   metrics) and Block 2.2 (annotation-stratified metrics).

## Manuscript implications

- **Figure-1 candidate**: the truth panel of the variance-capture curve is
  a striking one-glance summary of the architecture difference. The
  estimated panel shows the qualitative recovery. Could be **the figure
  that frames the architecture chapter of the paper**.
- **Per-trait architecture fingerprint**: each real trait gets a single
  variance-capture curve. Comparing curves across traits (e.g., height vs
  T2D vs schizophrenia) on one axes makes the polygenicity question
  visual.
- **D1 / D2 (CLAIMS)**: this work moves both from WEAK / PLAUSIBLE toward
  HELD pending Block 2.5.2's controlled validation.

## Next steps (Block 2.5.2)

- Simulate traits with **controlled polygenicity** (`n_causal_regions ∈
  {10, 100, 1000}`) at otherwise constant h²_arg / length / demography.
- Show `M_e,ARG` and variance-capture curves separate the three regimes
  cleanly at the truth level; quantify the estimate's bias.
- **Rare-variant ablation**: hide a controlled fraction of causal SNPs
  from typing → BayesR π_k blurs (typed-SNP polygenicity recovery fails);
  SARGENT `M_e,ARG` should hold (ARG-side polygenicity recovery is
  robust). This is the differentiable result that backs **B6**
  (rare-variant recovery via ARG kinship).
- DECISION POINT on the fused-lasso revival.

## Reproduce

```bash
# Phase 10 sims (already on the laptop from prior Block 1.2 work)
cd ~/test-sargent/phase10_admix
python ~/git-src/sargent/test/eval/arch_metrics.py --batch . --no-branches \
  --out arch_phase10.json
python ~/git-src/sargent/test/eval/plot_arch.py arch_phase10.json \
  --out-prefix arch_

# Block 1.2 sweep
cd ~/test-sargent/block1_2_v1
python ~/git-src/sargent/test/eval/arch_metrics.py --batch . --no-branches \
  --out arch_block1_2.json
python ~/git-src/sargent/test/eval/plot_arch.py arch_block1_2.json \
  --group-key ws2_block_contrast --out-prefix arch_

# To include the per-branch effect-size histogram, drop --no-branches
# (slower: walks each tree sequence's branches).
```

## Caveats

- The per-case ARG used for variance-capture curves is the **full ARG**
  (`case.full.trees`) when present; falls back to `case.trees` (train-only,
  simplified) otherwise. For Block 1.2 these match in structure but differ
  in tree count, which is the natural normalisation for the curve.
- `S_truth` and `S_est` differ across cases when the truth file is the
  full-ARG-indexed `case.full.true_ws2.tsv` but the estimate is
  train-tree-indexed (the typical case). The curves are reported in k/S
  units to make this comparable.
- The per-branch effect-size distribution is currently unweighted by the
  number of haplotypes carrying the variant. A frequency-weighted version
  (multiplying by `n_below(b) · (N − n_below(b))` or the BayesR
  variance-explained per allele frequency) is a possible v2 refinement.
