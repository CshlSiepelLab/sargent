# Cross-ancestry prediction under clumped per-segment ws² — comparing our result to Nurdan's

Author: Adam Siepel (+ Claude, laptop session)
Date: 2026-06-08

This report describes the SARGENT experiment we ran (Block 1.1 v3a) and
its result, so Nurdan can diff her setup against ours and we can find
where the two experiments disagree.

## TL;DR

We ran a cross-ancestry prediction sweep with a **clumped per-segment
ws² truth** (5% of trees active in a few contiguous blocks, 1000× higher
than the floor) and asked whether **per-segment ŵ² estimation** beats
**ignoring ws² entirely** (i.e. fitting a single uniform random
component). Across 30 cells spanning F_ST ∈ {0.07, 0.16, 0.30} and mix
fractions {0.00, 0.20}, the answer was **no**: the per-segment arm ties
the uniform arm (paired gap **+0.007 ± 0.035 r²_calib**, signs mixed
across cells).

Nurdan's preliminary work suggests per-segment estimation **does** help
cross-ancestry, so we want to identify where the two experiments differ.
This document lays out exactly what we did so she can check.

---

## The two arms we compared

We use SARGENT's CLI conventions; the equivalent MKLMM names are
listed for reference.

| arm | what it estimates | exact invocation |
|---|---|---|
| **Per-segment** (SARGENT-W) | per-tree ŵ²_s, fused-lasso TV smoothing on log ws², **no β fixed effect** | `sargent --arg-only --ws2-l2 0.5 ...` |
| **No-per-segment** (SARGENT-U / ARG-BLUP-uniform) | single uniform τ² across all trees (REML), no per-segment heterogeneity | computed analytically in `predict.py` from Σ_s G_s (or equivalently `sargent --arg-only --uniform-arg`) |

Both arms fit the **same model class** (y = u + ε with u ~ N(0, Σ_s w²_s G_s) and ε ~ N(0, σ²I)); only the prior/constraint on w²_s differs (per-tree-estimated vs. uniform).

MKLMM analogs (per `nurdan_compare/README.md`):
- Per-segment ≈ MKLMM **ARG-LMM** (AIREML per-kernel) or **Just-Random** (multi-kernel, no fixed effects).
- No-per-segment ≈ MKLMM **GRM-LMM (ARG)** (single averaged ARG kernel).

The cleanest 1:1 contrast on her side is therefore:
**ARG-LMM (or Just-Random) cross-ancestry r² vs. GRM-LMM(ARG) cross-ancestry r²**, on the same simulations.

---

## Our experiment — Block 1.1 v3a

### Demography

Simple 2-population split via `msprime` (driver `simulate.py
--demography transfer`). One ancestral pop, splits at a tunable time
into pop A (training) and pop B (testing). F_ST is controlled by the
split time:

| split (gen) | F_ST (Hudson) |
|---|---|
| 1500 | 0.07 |
| 3500 | 0.16 |
| 7000 | 0.30 |

Ne = 1e4 throughout. **No realistic demography** (no stdpopsim / Gravel
OOA / AmericanAdmixture in v3a — those were run separately in v3 with
the default *scattered* architecture and are NOT the experiment under
discussion here).

### Sample sizes and sequence

- `N_TRAIN = 500` (split between donor and target as `(1−mix)·N_TRAIN`
  from pop A and `mix·N_TRAIN` from pop B)
- `N_TEST = 200` (always from pop B)
- `LENGTH = 1_000_000` bp (1 Mb)
- `causal_fraction = 0.01`

### Heritabilities

- `h2_snp = 0.30` (BayesR-attributable typed-SNP variance)
- `h2_arg = 0.30` (ARG random-effect variance)
- Total broad-sense h² ≈ 0.60

The trait has both a SNP fixed-effect component (h²_snp) and an
ARG random-effect component (h²_arg). The per-segment arm
**ignores the SNP component** (because `--arg-only` sets β = 0)
and tries to absorb whatever it can via the ARG random effect.
**Important interaction**: with h²_snp > 0, omitting β is a
deliberate handicap on the per-segment arm; we found in v2 that
**including** β was worse cross-ancestry because β is the
cross-ancestry liability (the SNP effects don't port). The W arm
is the rescue test "does per-segment ŵ² help once the β liability
is removed?". If Nurdan's experiment has **h²_snp = 0** or uses a
SNP fixed effect that doesn't suffer cross-ancestry, that's one of
the dimensions to check.

### Clumped ws² architecture (the "clustered" truth)

Simulated via:
```
--ws2-architecture sparse
--ws2-n-blocks 3                # default; not set explicitly
--ws2-block-frac 0.05
--ws2-block-contrast 1000
```

Concretely (from `simulate.py:build_ws2_architecture`):
- The genome's ~S trees are arranged in their physical order along the
  sequence.
- `block_frac × n_trees` trees (5%) are flagged as "active"; the rest
  carry a floor ws² of `1 / contrast = 0.001`.
- The active trees are arranged into `n_blocks = 3` evenly-spaced
  **contiguous runs** along the sequence.
- The block trees carry ws² ≈ 1; floor trees carry ws² ≈ 0.001 → block-
  to-floor ratio = 1000×.
- Absolute scale is then renormalised so realised `var(u) = h²_arg`.

Causal variants are drawn per the standard `simulate.py` flow; the
realised per-tree contribution to y is u_s ~ N(0, ws²_s · h²_arg / Σ
ws²_s · G_s). In words: **the high-ws² blocks carry essentially all the
ARG variance, the floor trees carry almost none**.

### Typed-SNP coverage

`--bed-coverage 0.70` — 70% of variants in the BED panel.

(Held fixed across all 30 v3a cells. We do not vary typed coverage
within v3a; the v2/v3 sweeps that varied it kept the default scattered
architecture, so we don't have a 4D coverage × architecture sweep.)

### Cross-ancestry mix (target-population injection)

`MIX_FRAC_B ∈ {0.00, 0.20}`. mix = 0 means pure pop-A training; mix =
0.20 means 20% of training individuals come from pop B (the test pop).

### Replication

5 seeds per cell. Grid: 3 splits × 2 mix × 1 typed × 5 seeds = **30
cells**, 0 failures.

### Inference flags (the part most likely to differ vs. Nurdan)

Per-segment arm (W):
```
sargent --arg-only --ws2-l2 0.5 \
        --nprobe 10 --silent --outroot $dir/est_W \
        $dir/case.phen $dir/case.trees $dir/case
```
- `--arg-only`: drops β fixed effect (Stage-1 BayesR is skipped, β = 0).
- `--ws2-l2 0.5`: fused-lasso TV penalty on log ws², coefficient 0.5.
  This is **already the natural match to a contiguous-block
  architecture**: it penalises log-ws²[s+1] − log-ws²[s] which keeps
  adjacent trees coherent. **But** we did NOT use the sparsity
  penalty `--ws2-l1` which would have been the natural match to a
  5%-active mask. Setting `--ws2-l1 > 0` is the "v3b" test we have
  not run yet.

No-per-segment arm (U): no sargent call; the prediction is computed
analytically in `predict.py` as

> ŷ_test = (Σ_s G_s_test_train) · V_train⁻¹ · y_train,
> V_train = σ²I + τ² · (Σ_s G_s_train_train),

with τ² fit by REML on (y_train, K = Σ_s G_s).

### Metric

**r²_calib = squared Pearson correlation** between predicted and true
test phenotype. We learned in v2 that raw R² is destroyed cross-
ancestry by a single rare-in-train / common-in-test SNP that blows up
the BayesR-arm scale; r²_calib is calibration-invariant and is what
all our verdict tables report. (The driver also reports raw R²; the
two agree on within-pop arms and diverge dramatically on cross-pop
arms with a fixed-effect component.)

---

## Our result (Block 1.1 v3a)

Mean r²_calib over 5 seeds per cell; **per-segment** = SARGENT-W,
**no-per-segment** = SARGENT-U.

| F_ST | mix | no-per-segment (U) | per-segment (W) | W − U |
|---|---|---|---|---|
| 0.07 | 0.00 | 0.117 | 0.140 | +0.023 |
| 0.07 | 0.20 | 0.209 | 0.230 | +0.021 |
| 0.16 | 0.00 | 0.063 | 0.039 | −0.023 |
| 0.16 | 0.20 | 0.230 | 0.230 | +0.001 |
| 0.30 | 0.00 | 0.038 | 0.052 | +0.014 |
| 0.30 | 0.20 | 0.269 | 0.274 | +0.005 |

**Pooled (n = 30): W − U = +0.007 ± 0.035.**

Verdict: neutral. Per-segment ŵ² confers no positive cross-ancestry
benefit over uniform under this configuration. (For context: under
the *scattered* ws² truth tested in v2, per-segment was a clear drag
at −0.012 to −0.028 across F_ST. So clumping does remove the drag —
just not far enough to **beat** uniform.)

Per-cell raw data: `test/eval/block1_1_v3a_results.csv`.

---

## What we did NOT test

Honest list, because Nurdan's result may live in a regime we haven't
probed:

1. **Clumped architecture + realistic demography.** v3 (Gravel OOA,
   AmericanAdmixture, Ragsdale archaic) was run with the **default
   scattered** architecture, not the clumped one. We do not have a
   crossover where both stresses are present.
2. **`--ws2-l1` sparsity penalty.** The natural match to a 5%-active
   architecture is the sparsity-on-levels penalty (shrinks inactive
   trees to a floor), not the TV smoothing alone. We used `--ws2-l2
   0.5` only. A tuned (`--ws2-l1`, `--ws2-l2`) pair on the 30 existing
   cells is cheap and is on our queue as "v3b" — not yet run.
3. **Architecture parameter sweep.** We tested one cell of the
   (n_blocks × block_frac × contrast) space: (3, 0.05, 1000). We have
   not varied these. If Nurdan's clumping is finer (more blocks,
   smaller frac) or coarser (one big block), the per-segment-vs-
   uniform contrast could move.
4. **Higher per-block contrast.** 1000× was our pick; we have not
   tested 10× or 10000×.
5. **Pure h²_arg traits (h²_snp = 0).** Our trait has h²_snp = 0.3,
   h²_arg = 0.3. If Nurdan's trait is pure-ARG (no typed SNP signal),
   the W vs U comparison may behave differently because the per-
   segment arm doesn't need to ignore a β liability.
6. **Larger N_train.** 500 is modest for cross-ancestry. We have not
   tested whether per-segment ŵ² gains relative power at N = 2000 or
   5000.

---

## Dimensions to diff

The most likely places the two experiments differ. We're not
prescribing values; just listing the axes Nurdan can check her config
against ours.

| dimension | our v3a | Nurdan's? |
|---|---|---|
| Demography | msprime 2-pop split, Ne=1e4 | ? |
| F_ST (donor-target) | {0.07, 0.16, 0.30} | ? |
| N_train | 500 | ? |
| N_test | 200 | ? |
| Train mix fraction | {0.00, 0.20} | ? |
| Sequence length | 1 Mb | ? |
| Number of trees S | depends on sim (~few thousand) | ? |
| h²_snp | 0.30 | ? |
| h²_arg | 0.30 | ? |
| Causal fraction | 0.01 | ? |
| Typed-SNP coverage | 0.70 | ? |
| ws² architecture | 3 contiguous blocks, 5% active, 1000× contrast | ? (n_blocks? frac? contrast?) |
| ws² block placement | evenly spaced along the sequence | ? |
| Per-segment arm | `sargent --arg-only --ws2-l2 0.5` (TV-smoothed, no β, no l1) | likely MKLMM ARG-LMM with AIREML-per-kernel, no penalty |
| No-per-segment arm | analytic ARG-BLUP-uniform from Σ_s G_s (single τ²) | likely MKLMM GRM-LMM(ARG) — single averaged ARG kernel |
| Topology deduplication | none (all S trees used by SARGENT) | MKLMM deduplicates to unique topologies |
| Per-segment penalty | `--ws2-l2 0.5` (TV); no `--ws2-l1` | unpenalised AIREML |
| Per-segment parameterisation | variance ŵ²_s | MKLMM uses SD dr_s, V = σ²I + Σ dr²_s K_s |
| Metric | r²_calib (squared Pearson r) on test set | ? |
| Cross-ancestry direction | train pop A → test pop B (B-only test set) | ? |

The **two dimensions we'd most expect to matter** are:

(a) **Architecture parameters.** Our pick (3 blocks, 5% active, 1000×)
is one cell of a wide space. The signal-to-noise on per-segment ŵ² is
highly sensitive to whether each block has enough trees to support a
stable per-tree posterior. If Nurdan's clumps are larger (e.g., 1 block
spanning 20% of the genome at 50× contrast), the per-segment
estimator has more to work with per active tree, and the U baseline is
a worse approximation because the active region is concentrated.

(b) **Per-segment penalty / parameterisation.** Our W arm runs with
fused-lasso TV at λ=0.5 but no sparsity. MKLMM's AIREML on
deduplicated unique topologies is in a different regularisation
regime — it can drive per-kernel weights to zero individually but
doesn't enforce contiguity. With clumped truth and dedup, an unpenalised
AIREML may concentrate signal more sharply on the active topologies
than a fused-smoothed SARGENT does, recovering more of the clump-vs-
floor structure. That could be a real methodological difference
worth understanding.

---

## Concrete questions for Nurdan

1. **What architecture** did you simulate? Specifically:
   - How many active blocks?
   - What fraction of trees / sequence length is active?
   - What is the block-to-floor ws² contrast?
   - Are the blocks contiguous along the sequence?
2. **What demography**? msprime 2-pop split, msprime admixture,
   stdpopsim, real haplotype reference, …?
3. **What N_train, N_test, and what cross-ancestry direction** (train
   on which pop, test on which)?
4. **What sequence length** and **how many trees** does the per-
   segment estimator see?
5. **What metric** are you reporting? Raw R² on a held-out cohort, OOS
   prediction r, calibration-invariant squared correlation, h² recovery
   on the test cohort, …?
6. **For the no-per-segment baseline**: is it GRM-LMM (single averaged
   ARG kernel) or something else (typed-SNP BLUP, BayesR-only, etc.)?
   We want the same baseline.
7. **Is the trait pure-ARG (h²_snp = 0) or does it include typed-SNP
   variance**? Our v3a has h²_snp = 0.3 alongside h²_arg = 0.3, and
   our W arm drops β to avoid the cross-ancestry β liability — if her
   trait is pure-ARG, dropping β is a non-issue and the per-segment
   arm has nothing to "rescue".
8. Can you share one simulation case directory (or its `params.json`)
   from the experiment that shows the per-segment win? With one
   concrete case we can run **our** estimator on **your** sim and
   isolate whether the gap is in the data or in the estimator.

---

## Reproduce our v3a

```bash
# From the sargent repo root, with the sargent mamba env active.
BASE=$HOME/test-sargent/block1_1_v3a_repro
SIM_EXTRA='--ws2-architecture sparse --ws2-block-frac 0.05 --ws2-block-contrast 1000' \
SPLITS='1500 3500 7000' \
MIX_FRAC_B='0.00 0.20' \
TYPED_FRACS='0.70' \
SEEDS='1 2 3 4 5' \
N_TRAIN=500 N_TEST=200 LENGTH=1000000 \
H2_SNP=0.3 H2_ARG=0.3 CAUSAL_FRAC=0.01 \
BASE=$BASE bash test/eval/portability.sh

# The W arm (per-segment, --arg-only --ws2-l2 0.5) is added by:
BASE=$BASE VARIANTS=W bash test/eval/sargent_variants_v2.sh

# Aggregate:
python test/eval/plot_portability.py $BASE
```

Per-cell results land at `$BASE/sp<split>_mix<m>_typed<t>_s<seed>/{est,est_W}.predict.json`; aggregate at `$BASE/portability_table.md`.

---

## References (internal)

- v3a design + verdict: `test/eval/block1_1_design.md` § "Block 1.1 v3a — Clustered-architecture cross-ancestry rescue (2026-06-05)" (lines 638-699).
- v3a per-cell CSV: `test/eval/block1_1_v3a_results.csv`.
- Driver: `test/eval/portability.sh` (simulation + main fit + predict); `test/eval/sargent_variants_v2.sh` (adds the W arm).
- Simulator: `scripts/simulate.py:build_ws2_architecture` for the clumped truth.
- Predictor: `test/eval/predict.py` (4-arm OOS prediction including ARG-BLUP-uniform analytic).
- Nurdan-side overview: `test/eval/nurdan_compare/nurdan-overview.md` + `README.md`.
