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

- [ ] `test/eval/localization_sweep.sh`: driver + aggregator. Reuses
      `oracle_ws2.py` as the per-cell runner (it already computes everything
      we want).
- [ ] Smoke-test on a 2-cell run.
- [ ] Hold v1 execution behind Adam's go.
