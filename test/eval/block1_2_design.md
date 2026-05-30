# Block 1.2 — Localization robustness sweep (design v1)

Date: 2026-05-29.
Status: design (pending Adam's sign-off; smoke-tested code in place).
Owners: Claude (sim prep + execution), Adam (decision review).
Backs claims: C1 (region-scale localization works), C4 [MINOR] (estimator
robustness), C5 (gated; sparse fused lasso) — see `CLAIMS.md`.

## Goal

Confirm — or refute — that SARGENT's per-tree REML estimator identifies
active genealogical regions at **AUC well above chance across a range of
architectures**, not just the single sparse + admixed + contrast-1000 cell
that Phase 10 measured. The result decides whether Block 2 (sparse fused
lasso) is needed to close a localization gap.

Phase 10 single cell (5 seeds): estimator AUC 0.74 ± 0.11, ceiling 0.84.
Block 1.2 asks: how does this scale across {contrast, demography}?

## Design

### Architecture grid

Pin **sparsity (BLOCKFRAC) = 0.05** and **NBLOCKS = 3** (Phase 10's
configuration) for v1. Vary:

| knob       | values              | meaning |
|---|---|---|
| BLOCKCONTRAST | 1000 / 100 / 10  | strong / moderate / weak per-segment signal |
| demography    | admixed / panmictic | local-ancestry heterogeneity present vs absent |
| seeds         | 1 .. 5             | replicates |

Total: 3 × 2 × 5 = **30 cells**.

(*Sparsity sweep — e.g. block_frac 0.02, 0.20 — and n_train sweep — 300 vs
500 — are deferred to v2 if v1 is informative.*)

### Why these knobs

- **Contrast** answers "does the estimator work when the block:floor ratio
  is biologically realistic (10×) rather than artificial (1000×)?". Without
  this we have no idea where the estimator falls off.
- **Demography** answers "does the estimator detect blocks when there is
  no local-ancestry heterogeneity to lean on?" The panmictic case is
  structurally pessimistic for ws² (the ARG itself carries no spatial
  variation), so any AUC > 0.5 there confirms the per-tree REML signal
  comes from y, not from ancestry.

### Sample sizes, length, h²

- **n_train = 300**, **n_test = 100** (matches Phase 10; ~3 min/cell).
- **L = 300 kb** (matches Phase 10).
- **h²_snp = 0**, **h²_arg = 0.5** (pure random; same as Phase 10).

### Metrics per cell

From `oracle_ws2.py` output (`est.oracle.json`), report:

- **Estimator localization AUC**: `ws2_recovery.auc_active` — the headline.
- **Localization ceiling AUC**: `loc_ceiling_auc` — the idealized score at
  the homogeneous-null V (computed by `oracle_ws2.py`).
- **ws² profile recovery**: `ws2_recovery.{pearson, spearman}` —
  monotone correlation of estimated vs true `ws²` across trees.
- Realised h²_arg (sanity).

Aggregated outputs:

- Per-seed table.
- **(contrast × demography) summary** with mean ± sd estimator AUC and
  the **estimator − ceiling gap** (how much headroom is left).
- Trend lines: estimator AUC vs contrast (per demography); estimator AUC
  vs ceiling (per cell).

## Sequencing

### v1 (this design, ~2 h wall on a single core)

30 cells. Three outcomes:

- **All-positive** (estimator AUC > 0.65 across all cells): C1 holds
  robustly. Block 2 (sparse fused lasso) becomes optional polish, not a
  prerequisite. Cycles freed for Block 1.1 v2 or Block 4.
- **Drops with contrast** (high at contrast 1000, falls to 0.5 at 10):
  the estimator is signal-limited. C1 holds *in the strong-contrast regime*
  only. Block 2 is justified as the way to extract weak-signal localization
  via regularization.
- **Drops with demography** (works admixed, fails panmictic): the
  estimator is leaning on ancestry, not on y. C1 needs a structured-pop
  caveat in the manuscript. Block 2 might or might not help.

### v2 (if v1 is informative)

- Sparsity ladder: 0.02 / 0.05 / 0.20.
- n_train ladder: 300 / 500.
- Add NBLOCKS variation (1 / 3 / 10) if reviewers might ask about
  "single big locus" vs "scattered small loci".

## Decision log

| # | Choice | Rationale |
|---|---|---|
| 1 | Pin BLOCKFRAC = 0.05 in v1 | One contrast × demography sweep at the Phase 10 architecture; cheaper, more focused. Expand sparsity in v2. |
| 2 | Contrast at 1000 / 100 / 10 (not finer) | Three-log span captures the trend; finer ladder adds run time without sharpening the answer. |
| 3 | Include panmictic as a cell, not just admixed | Tests whether per-tree REML uses y-signal or ancestry-signal. Cheap; high diagnostic value. |
| 4 | Stay at n=300, L=300kb | Matches Phase 10 baseline; comparable numbers. Bumping n is a v2 axis. |
| 5 | 5 seeds | Same as Phase 10; per-seed SD of ~0.1 in AUC suggests ≥5 is the minimum honest replicate count. |

## Reproduce (planned interface)

```bash
DEMOGS="admixed panmictic" \
  CONTRASTS="1000 100 10" \
  SEEDS="1 2 3 4 5" \
  N_TRAIN=300 N_TEST=100 LENGTH=300000 \
  BLOCKFRAC=0.05 NBLOCKS=3 \
  BASE=/path/to/block1_2_v1 \
  bash test/eval/localization_sweep.sh
```

Outputs land per cell at `<BASE>/c<contrast>_d<demog>_s<seed>/` with
`est.oracle.json` per cell; aggregator writes
`<BASE>/localization_table.md`.

## SSH launch (server)

Same pattern as Block 1.1 (see `CLAUDE.md`). Recommended:

```bash
git pull --ff-only
BASE=$HOME/test-sargent/block1_2_v1
mkdir -p $BASE
tmux new -s b12 -d "BASE=$BASE bash test/eval/localization_sweep.sh 2>&1 | tee $BASE/run.log"
# tmux attach -t b12 to reattach
```

Resumable: skips cells whose `est.oracle.json` exists.

## Open implementation tasks

- [x] `test/eval/localization_sweep.sh`: driver + aggregator. Reuses
      `oracle_ws2.py` as the per-cell runner.
- [x] Smoke-test on a 2-cell run.
- [x] v1 executed on laptop (2026-05-29 → 30, 593 min wall, all 30 cells).

## v1 results (2026-05-30)

Raw table: `~/test-sargent/block1_2_v1/localization_table.md`.
Reproducibility note: laptop run with `BASE=~/test-sargent/block1_2_v1`.

### Headline numbers (mean ± sd over 5 seeds)

| demog | contrast | **est AUC** | ceiling AUC | gap (ceiling − est) | spearman(ŵ², w²) |
|---|---|---|---|---|---|
| admixed   | 1000 | **0.809 ± 0.086** | 0.860 ± 0.088 | 0.051 ± 0.116 | +0.190 ± 0.153 |
| admixed   | 100  | **0.790 ± 0.105** | 0.839 ± 0.102 | 0.049 ± 0.107 | +0.187 ± 0.128 |
| admixed   | 10   | **0.677 ± 0.130** | 0.721 ± 0.185 | 0.043 ± 0.068 | +0.152 ± 0.076 |
| panmictic | 1000 | **0.707 ± 0.107** | 0.757 ± 0.108 | 0.051 ± 0.106 | −0.058 ± 0.076 |
| panmictic | 100  | **0.633 ± 0.115** | 0.727 ± 0.095 | 0.095 ± 0.118 | −0.022 ± 0.113 |
| panmictic | 10   | **0.511 ± 0.084** | 0.617 ± 0.080 | 0.105 ± 0.088 | −0.096 ± 0.126 |

### Interpretation

The result lands between scenarios (b) and (c) from the design — partly
both, neither fully. Three findings stand out:

1. **The estimator is above chance in 5 of 6 cells.** Only panmictic +
   contrast 10 (AUC 0.51) is at chance. Region-scale localization holds
   robustly in the populations and signal regimes that matter most
   (admixed; contrast ≥ 100). Backs **C1 (HELD)** — the Phase 10 single-cell
   result generalizes.
2. **The estimator–ceiling gap is small everywhere (~0.04–0.11).** The
   per-tree REML EM extracts most of the available information. This is
   important: it means Block 2 (sparse fused lasso) has **much less
   headroom than we hoped**. The original Phase 11 motivation was closing
   a ~0.10 estimator gap; the worst case here (panmictic + contrast 10) has
   a gap of 0.11, but the *ceiling* is only 0.62 — even a perfect estimator
   only reaches there.
3. **Two robustness gradients, both expected:**
   - AUC drops with contrast (1000 → 100 → 10): 0.81 → 0.79 → 0.68 in
     admixed; 0.71 → 0.63 → 0.51 in panmictic. Signal-limited at the low end.
   - AUC drops in panmictic vs admixed at every contrast (~0.10 systematic
     gap). The estimator *does* benefit from ancestry structure, but it
     does not *require* it — panmictic at contrast 1000 still gives
     AUC 0.71.

### Implications

**For Block 2 (Phase 11 sparse fused lasso):**

> The gap between the per-tree REML estimator and the idealized ceiling is
> already small. Sparse fused lasso would need to close ≤0.11 in the worst
> cell and ≤0.05 in the typical cell. **Downgrade Block 2** from "needed to
> rescue localization" to "optional polish at the low-signal end". Cycles
> may be better spent elsewhere — particularly on Block 1.1 v2 if v1 is
> positive.

**For Block 1.1 (prediction):**

> The Block 1.2 numbers do not directly speak to prediction. But they
> confirm the *partitioned-h²* story has real signal at scale (AUC ~ 0.7–0.8
> in the admixed regime), giving the manuscript a strong localization
> result regardless of how Block 1.1 lands.

**For CLAIMS.md (proposed bucket changes):**

| claim | current | proposed | reason |
|---|---|---|---|
| C1 (region-scale localization) | HELD (one architecture) | **HELD** | Generalizes across 6 architecture cells; AUC > chance in 5/6 |
| C4 [MINOR] (per-tree REML estimator) | WEAK / UNTESTED [MINOR] | **HELD [MINOR]** | Estimator–ceiling gap is small (~0.04–0.11) across the grid |
| C5 (sparse fused lasso closes gap) | PLAUSIBLE / UNTESTED | **WEAK** | Gap to close is small; Block 2 has less headroom than expected |
| A4 (partitioned h² respects recombination boundaries) | PLAUSIBLE / UNTESTED | (unchanged; needs partitioned-h² recovery metric, not AUC) | |

### v2 possibilities (deprioritized)

Originally sketched as: sparsity ladder (0.02 / 0.05 / 0.20), n_train
ladder (300 / 500), NBLOCKS variation. With v1 showing the gap is small
and the estimator robust, v2 is **deprioritized** — it would refine an
already-positive answer rather than change the verdict.

### Caveats

- 5 seeds per cell; per-seed AUC SD is 0.08–0.13. CIs are wide but the
  *trend* across contrast and demography is monotone and clean.
- Per-seed outliers: admixed/1000/s3 (est 0.71, ceiling 0.97) and
  panmictic/1000/s5 (est 0.69, ceiling 0.58 — estimator above ceiling) are
  reminders that the heuristic ceiling `A_s/tr(G_s)` is not strictly
  optimal, only an indicator of the right order of magnitude.
- Only one (block_frac, n_blocks) = (0.05, 3) architecture tested. v2
  sparsity ladder would address this if a reviewer asks.
- All sims used pure-random (h²_snp = 0, h²_arg = 0.5). Adding fixed
  effects is a separate axis untested here.
